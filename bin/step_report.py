#!/usr/bin/env python3
"""Render one step PDF: page 1 plain language, page 2+ technical.

Usage: step_report.py SUMMARY.json STEP_TEXT.yaml OUT.pdf [--plots a.png b.png ...]

The summary JSON (schema mtb-step-summary/1) is written by each step:
  sample_id, run_id, step, status (PASS|WARN|FAIL), pipeline_version, generated,
  values   {name: value}                      -> fills {placeholders} in the wording file
  metrics  [{name, plain, value, unit, threshold, status, provisional}]
  tools    [{name, version, container, command}]
  files    {inputs: [{path, checksum}], outputs: [...]}
  messages [str]
  results  {title, columns, rows, badge_col}   optional result table on page 1 (badge words in RESULT)
  tables   [{title, columns, rows, widths}]    optional extra tables on the technical pages
"""
import argparse, json
import yaml
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate,
                                Spacer, Table, TableStyle)

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
NAVY = colors.HexColor("#1F3A5F")
# status colours are never the only signal: the word PASS / CHECK / FAIL is always printed
STATUS = {"PASS": ("PASS", "#0a7d32"), "WARN": ("CHECK", "#b35c00"), "FAIL": ("FAIL", "#b3261e"),
          "INFO": ("INFO", "#52514e")}
# result badges (findings, not QC states); the word is always printed
RESULT = {"DETECTED": "#b3261e", "NOT DETECTED": "#0a7d32", "INDETERMINATE": "#52514e",
          "MINORITY": "#b35c00", "REVIEW": "#b35c00"}

ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=ss["Title"], fontSize=17, leading=21, alignment=0, spaceAfter=2,
                    textColor=colors.HexColor(INK))
SUB = ParagraphStyle("SUB", parent=ss["Normal"], fontSize=8.5, textColor=colors.HexColor(INK2))
H2 = ParagraphStyle("H2", parent=ss["Heading2"], fontSize=12, textColor=NAVY, spaceBefore=8, spaceAfter=3)
BIG = ParagraphStyle("BIG", parent=ss["Normal"], fontSize=11, leading=15, spaceAfter=6)
B = ParagraphStyle("B", parent=ss["Normal"], fontSize=9.2, leading=12.4, spaceAfter=4)
C = ParagraphStyle("C", parent=B, fontSize=8, leading=10, spaceAfter=0)
CH = ParagraphStyle("CH", parent=C, textColor=colors.white, fontName="Helvetica-Bold")
MONO = ParagraphStyle("MONO", parent=C, fontName="Courier", fontSize=7.2, leading=8.8)


class NumberedCanvas(canvas.Canvas):
    """Adds 'page x of y' and the footer once the page count is known."""
    footer = ""

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self._pages = []

    def showPage(self):
        self._pages.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        n = len(self._pages)
        for st in self._pages:
            self.__dict__.update(st)
            self.setFont("Helvetica", 7)
            self.setFillColor(colors.HexColor(MUTED))
            self.drawString(15 * mm, 9 * mm, self.footer)
            self.drawRightString(195 * mm, 9 * mm, f"Page {self._pageNumber} of {n}")
            super().showPage()
        super().save()


def esc(x, limit=None):
    x = str(x)
    if limit and len(x) > limit:
        x = x[:limit] + " ... (truncated; full list in the step's data/ folder)"
    return x.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def badge(status, big=False):
    word, col = STATUS.get(status) or (status, RESULT.get(status.split(" (")[0], "#52514e"))
    st = ParagraphStyle("bd", parent=C, fontName="Helvetica-Bold", textColor=colors.white, alignment=1,
                        fontSize=14 if big else 7.5, leading=17 if big else 9)
    w = 34 * mm if big else max(15 * mm, len(word) * 1.9 * mm)
    t = Table([[Paragraph(word, st)]], colWidths=[w])
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(col)),
                           ("TOPPADDING", (0, 0), (-1, -1), 4 if big else 1),
                           ("BOTTOMPADDING", (0, 0), (-1, -1), 5 if big else 2)]))
    return t


def table(rows, widths, status_col=None):
    data = []
    for r, row in enumerate(rows):
        out = []
        for c, v in enumerate(row):
            if r and c == status_col:
                out.append(badge(v))
            else:
                out.append(Paragraph(v if r == 0 else esc(v, 700), CH if r == 0 else C))
        data.append(out)
    t = Table(data, colWidths=widths, repeatRows=1)
    st = [("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#c3c2b7")), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
          ("BACKGROUND", (0, 0), (-1, 0), NAVY), ("LEFTPADDING", (0, 0), (-1, -1), 4),
          ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]
    t.setStyle(TableStyle(st))
    return t


def fmt(v, unit=""):
    if isinstance(v, float):
        v = f"{v:,.2f}"
    elif isinstance(v, int):
        v = f"{v:,}"
    return f"{v} {unit}".strip()


def header(s, title, summ):
    s.append(Paragraph(esc(title), H1))
    s.append(Paragraph(f"Sample <b>{esc(summ['sample_id'])}</b> &nbsp;|&nbsp; Run {esc(summ['run_id'])} "
                       f"&nbsp;|&nbsp; Step {esc(summ['step'])} &nbsp;|&nbsp; {esc(summ['generated'])} "
                       f"&nbsp;|&nbsp; Pipeline {esc(summ['pipeline_version'])}", SUB))
    s.append(Spacer(1, 6))


def build(summ, text, out, plots):
    vals = {**summ.get("values", {}), "sample_id": summ["sample_id"]}
    status = summ["status"]
    s = []

    # ---------------- page 1: plain language
    header(s, text["title_plain"], summ)
    s.append(Table([[badge(status, big=True),
                     Paragraph(esc(text["result"][status].strip().format(**vals)), BIG)]],
                   colWidths=[40 * mm, 140 * mm],
                   style=[("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (0, 0), 0)]))
    s.append(Paragraph("What this step does", H2))
    s.append(Paragraph(esc(text["what"].strip()), B))
    s.append(Paragraph("Why it matters", H2))
    s.append(Paragraph(esc(text["why"].strip()), B))
    res = summ.get("results")
    if res:
        s.append(Paragraph(esc(res["title"]), H2))
        n = len(res["columns"])
        widths = res.get("widths_mm") or [180 / n] * n
        s.append(table([res["columns"]] + res["rows"], [w * mm for w in widths], status_col=res.get("badge_col")))
    s.append(Paragraph("What we checked", H2))
    rows = [["What we checked", "Result for this sample", "OK?"]]
    # reference-only (INFO) numbers stay on the technical pages
    rows += [[m["plain"], fmt(m["value"], m.get("unit", "")), m["status"]] for m in summ["metrics"]
             if m["status"] != "INFO"]
    s.append(table(rows, [80 * mm, 75 * mm, 25 * mm], status_col=2))
    s.append(Paragraph("What happens next", H2))
    s.append(Paragraph(esc(text["next"][status].strip()), B))
    s.append(PageBreak())

    # ---------------- page 2+: technical
    header(s, text["title_tech"], summ)
    s.append(Paragraph("Method", H2))
    s.append(Paragraph(esc(text["technical_method"].strip().format(**vals)), B))
    s.append(Paragraph("Metrics against thresholds", H2))
    rows = [["Metric", "Value", "Threshold (conf/thresholds.yaml)", "Status"]]
    for m in summ["metrics"]:
        th = m.get("threshold", "-") + (" (provisional)" if m.get("provisional") else "")
        rows.append([m["name"], fmt(m["value"], m.get("unit", "")), th, m["status"]])
    s.append(table(rows, [55 * mm, 40 * mm, 65 * mm, 20 * mm], status_col=3))
    for p in plots:
        s.append(Spacer(1, 6))
        iw, ih = ImageReader(p).getSize()
        s.append(Image(p, width=180 * mm, height=180 * mm * ih / iw))
    for tb in summ.get("tables", []):
        s.append(Paragraph(esc(tb["title"]), H2))
        n = len(tb["columns"])
        widths = tb.get("widths_mm") or [180 / n] * n
        s.append(table([tb["columns"]] + tb["rows"], [w * mm for w in widths] if tb["rows"] else None)
                 if tb["rows"] else Paragraph("None.", B))
    if summ.get("messages"):
        s.append(Paragraph("Messages", H2))
        for m in summ["messages"]:
            s.append(Paragraph("&bull; " + esc(m), B))
    s.append(Paragraph("Software", H2))
    rows = [["Tool", "Version", "Container"]] + [[t["name"], t["version"], t.get("container", "-")]
                                                 for t in summ["tools"]]
    s.append(table(rows, [35 * mm, 30 * mm, 115 * mm]))
    for t in summ["tools"]:
        if t.get("command"):
            s.append(Spacer(1, 3))
            s.append(KeepTogether([Paragraph(f"<b>Command ({esc(t['name'])})</b>", C),
                                   Paragraph(esc(t["command"]), MONO)]))
    s.append(Paragraph("Files and checksums", H2))
    rows = [["Role", "File", "Checksum"]]
    for role in ("inputs", "outputs"):
        rows += [[role[:-1], f["path"], f.get("checksum", "-")] for f in summ.get("files", {}).get(role, [])]
    s.append(table(rows, [16 * mm, 74 * mm, 90 * mm]))

    NumberedCanvas.footer = (f"Internal step report - not for clinical release  |  {summ['sample_id']}  |  "
                             f"step {summ['step']}")
    doc = SimpleDocTemplate(out, pagesize=A4, leftMargin=15 * mm, rightMargin=15 * mm, topMargin=14 * mm,
                            bottomMargin=16 * mm, title=f"{summ['sample_id']} {summ['step']}",
                            author=f"MTB_PIPELINE {summ['pipeline_version']}")
    doc.build(s, canvasmaker=NumberedCanvas)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("summary"); ap.add_argument("text"); ap.add_argument("out")
    ap.add_argument("--plots", nargs="*", default=[])
    a = ap.parse_args()
    build(json.load(open(a.summary)), yaml.safe_load(open(a.text)), a.out, a.plots)
