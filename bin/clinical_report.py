#!/usr/bin/env python3
"""Step 12: the clinical report (antibiogram + resistance details, per the project's approved
report design) and the merged all-steps PDF. Everything printed comes from the step summary JSONs and
the step 09 tables, so the report cannot disagree with the step PDFs.

Usage: clinical_report.py --summaries DIR --drug_results TSV --evidence TSV --variants TSV --regions BED
                          --drugs TSV --meta JSON --step_pdfs PDF... --lab_name S --operator S
                          --banner S --version V --out_dir .
"""
import argparse, csv, datetime, glob, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pypdf import PdfWriter
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table,
                                TableStyle)
from step_report import NumberedCanvas, esc
from mtbstep import INK, INK2

NAVY = colors.HexColor("#1F3A5F")
RES = {"DETECTED": "#b3261e", "NOT DETECTED": "#0a7d32", "INDETERMINATE": "#6b6a66"}
LETTER = {"DETECTED": "R", "NOT DETECTED": "S", "INDETERMINATE": "?"}
LIMS = [("client_id", "Client ID"), ("test_id", "Test ID"), ("patient_name", "Patient name"), ("age", "Age"),
        ("sex", "Sex"), ("referring_doctor", "Referring doctor"), ("sample_type", "Sample type"),
        ("sample_source", "Sample source"), ("receipt_date", "Sample receipt date"),
        ("sample_temperature", "Temperature on receipt"), ("sample_qc_date", "Sample QC date"),
        ("dna_260_280", "DNA 260/280"), ("dna_conc_ng_ul", "DNA concentration (ng/µl)")]
LIB = [("library_type", "Library type"), ("library_date", "Library date"), ("library_qc", "Library QC"),
       ("sequencer", "Genome sequencer"), ("run_date", "Run date")]

ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=ss["Title"], fontSize=15, leading=19, alignment=0, spaceAfter=1,
                    textColor=colors.HexColor(INK))
SUB = ParagraphStyle("SUB", parent=ss["Normal"], fontSize=8.3, textColor=colors.HexColor(INK2))
H2 = ParagraphStyle("H2", parent=ss["Heading2"], fontSize=11, textColor=NAVY, spaceBefore=7, spaceAfter=3)
B = ParagraphStyle("B", parent=ss["Normal"], fontSize=8.6, leading=11.2, spaceAfter=3)
C = ParagraphStyle("C", parent=B, fontSize=7.8, leading=9.6, spaceAfter=0)
CB = ParagraphStyle("CB", parent=C, fontName="Helvetica-Bold")
CH = ParagraphStyle("CH", parent=C, textColor=colors.white, fontName="Helvetica-Bold")
BANNER = ParagraphStyle("BAN", parent=B, fontName="Helvetica-Bold", textColor=colors.white, alignment=1,
                        backColor=colors.HexColor("#b3261e"), borderPadding=4, spaceAfter=6)
BIG = ParagraphStyle("BIG", parent=B, fontSize=12, leading=15, fontName="Helvetica-Bold")


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def kv_table(pairs, cols=2):
    """Label/value grid, `cols` pairs per row."""
    rows, row = [], []
    for k, v in pairs:
        row += [Paragraph(esc(k), CB), Paragraph(esc(v if v not in (None, "") else "—"), C)]
        if len(row) == 2 * cols:
            rows.append(row)
            row = []
    if row:
        rows.append(row + [""] * (2 * cols - len(row)))
    w = 180 / cols
    t = Table(rows, colWidths=[w * 0.42 * mm, w * 0.58 * mm] * cols)
    t.setStyle(TableStyle([("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#c3c2b7")),
                           ("BACKGROUND", (0, 0), (-1, -1), colors.white), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                           ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)] +
                          [("BACKGROUND", (c, 0), (c, -1), colors.HexColor("#F2F5F8")) for c in range(0, 2 * cols, 2)]))
    return t


def grid(rows, widths, badge_col=None):
    data = []
    for r, row in enumerate(rows):
        out = []
        for c, v in enumerate(row):
            if r and c == badge_col and v:
                st = ParagraphStyle("b", parent=CB, textColor=colors.white, alignment=1)
                b = Table([[Paragraph(esc(v), st)]], colWidths=[widths[c] * mm - 6])
                b.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(RES.get(v, "#6b6a66"))),
                                       ("TOPPADDING", (0, 0), (-1, -1), 1), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]))
                out.append(b)
            else:
                out.append(Paragraph(v if (r == 0 or v.startswith("<b>")) else esc(v), CH if r == 0 else C))
        data.append(out)
    t = Table(data, colWidths=[w * mm for w in widths], repeatRows=1)
    t.setStyle(TableStyle([("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#c3c2b7")),
                           ("BACKGROUND", (0, 0), (-1, 0), NAVY), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                           ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]))
    return t


def antibiogram(results, out):
    n = len(results)
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    cols = [RES[r["call"]] for r in results]
    wedges, _ = ax.pie([1] * n, colors=cols, startangle=90, counterclock=False,
                       wedgeprops={"width": 0.36, "edgecolor": "white", "linewidth": 2})
    import math
    for w, r in zip(wedges, results):
        ang = math.radians((w.theta1 + w.theta2) / 2)
        ax.text(0.82 * math.cos(ang), 0.82 * math.sin(ang), LETTER[r["call"]], ha="center", va="center",
                color="white", fontsize=9, fontweight="bold")
        ax.text(1.17 * math.cos(ang), 1.17 * math.sin(ang), r["code"], ha="center", va="center", color=INK,
                fontsize=8)
    ax.text(0, 0.08, "R  resistance mutation detected", ha="center", fontsize=6.8, color=INK2)
    ax.text(0, -0.06, "S  not detected", ha="center", fontsize=6.8, color=INK2)
    ax.text(0, -0.20, "?  indeterminate", ha="center", fontsize=6.8, color=INK2)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(out, dpi=220, transparent=False, facecolor="white")


def main():
    ap = argparse.ArgumentParser()
    for k in ("summaries", "drug_results", "evidence", "variants", "regions", "drugs", "meta", "lab_name",
              "operator", "banner", "version", "out_dir"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--step_pdfs", nargs="*", default=[])
    a = ap.parse_args()
    S = {os.path.basename(p).split(".")[0]: json.load(open(p))
         for p in glob.glob(os.path.join(a.summaries, "*.summary.json"))}
    meta = json.load(open(a.meta))
    sid, run = meta["id"], meta["run_id"]
    lims = meta.get("lims", {})
    results = tsv(a.drug_results)
    evidence = [e for e in tsv(a.evidence) if e["status"] == "used"]
    variants = {(v["pos"]): v for v in tsv(a.variants)}
    locus = {}
    for line in open(a.regions):
        f = line.rstrip("\n").split("\t")
        locus[f[3]] = f[4]

    def metric(step, name):
        return next((m for m in S[step]["metrics"] if m["name"] == name), {"value": "—", "status": "INFO"})

    # ---------- QC gates
    gates = [("Input files intact", "00_input", None),
             ("Read yield", "01_read_qc", "Raw yield"), ("Base quality (Q30)", "01_read_qc", "Q30 bases after trimming"),
             ("MTB complex (Mykrobe)", "02_species", "Mykrobe: MTB complex"),
             ("MTB complex (Sylph)", "02_species", "Sylph: MTB complex"),
             ("Reads mapped to H37Rv", "03_alignment", "Reads mapped to H37Rv"),
             ("Mean depth", "04_coverage", "Mean depth"), ("Genome at >= 10x", "04_coverage", "Genome at >= 10x")]
    gate_rows, qc_fail, qc_warn = [], False, False
    for label, step, name in gates:
        if name is None:
            st, val = S[step]["status"], "all checks"
        else:
            m = metric(step, name)
            st, val = m["status"], f"{m['value']} {m.get('unit', '')}".strip()
        qc_fail |= st == "FAIL"
        qc_warn |= st == "WARN"
        gate_rows.append([label, val, {"PASS": "Pass", "WARN": "Pass (see note)", "FAIL": "Fail"}.get(st, st)])
    qc = "FAIL" if qc_fail else ("PASS with comments" if qc_warn else "PASS")

    mtb = metric("02_species", "Mykrobe: MTB complex")["status"] == "PASS" and \
        metric("02_species", "Sylph: MTB complex")["status"] == "PASS"
    v02, v04, v08 = S["02_species"]["values"], S["04_coverage"]["values"], S["08_lineage"]["values"]
    v10, v11 = S["10_mixed_infection"]["values"], S["11_classification"]["values"]
    hdr = (S["00_input"].get("data") or {}).get("header") or {}

    s = []
    if a.banner:
        s.append(Paragraph(esc(a.banner), BANNER))
    s.append(Paragraph(esc(a.lab_name), SUB))
    s.append(Paragraph("Mycobacterium tuberculosis Whole-Genome Sequencing: Drug-Resistance Report", H1))
    s.append(Paragraph(f"Sample <b>{esc(sid)}</b> &nbsp;|&nbsp; Run {esc(run)} &nbsp;|&nbsp; Report date "
                       f"{datetime.date.today():%d-%m-%Y} &nbsp;|&nbsp; MTB_PIPELINE {esc(a.version)}", SUB))

    s.append(Paragraph("1. Sample summary", H2))
    s.append(kv_table([("Sample ID", sid)] + [(lbl, lims.get(k, "")) for k, lbl in LIMS]))
    s.append(Paragraph("2. Genome sequencing summary", H2))
    s.append(kv_table([(lbl, lims.get(k, "")) for k, lbl in LIB] +
                      [("Instrument ID (from reads)", hdr.get("instrument", "not in read header")),
                       ("Flow cell / lane", f"{hdr.get('flowcell', '—')} / {hdr.get('lane', '—')}" if hdr else "—"),
                       ("Sequence yield", f"{S['01_read_qc']['values'].get('raw_yield_gb', '—')} Gb")]))
    s.append(Paragraph("3. Genomic analysis summary", H2))
    s.append(kv_table([("Application", f"MTB_PIPELINE {a.version}"), ("Catalogue", "WHO mutation catalogue, 2nd ed. (2023)"),
                       ("Analysis date", f"{datetime.date.today():%d-%m-%Y}"),
                       ("Analysis run by", a.operator), ("QC check", qc),
                       ("Thresholds version", S["04_coverage"].get("thresholds_version", "—"))]))
    s.append(Paragraph("4. Clinical summary", H2))
    s.append(kv_table([("M. tuberculosis complex detected", "Yes" if mtb else "No / not confirmed"),
                       ("Species", v02.get("species", "—")),
                       ("Lineage", f"{v08.get('lineage', '—')} ({v08.get('family', '—')}), {v08.get('sub_lineage', '')}"),
                       ("Spoligotype", f"{v08.get('spoligo', '—')} (SIT {v08.get('sit', '—')})"),
                       ("Mean depth of coverage", f"{v04.get('mean_depth', '—')}x"),
                       ("Genome covered at >= 10x", f"{v04.get('breadth10', '—')}%")]))

    if qc_fail:
        s.append(Paragraph("5. Result", H2))
        s.append(Paragraph("QC FAILED: drug-resistance results are not reported for this sample. "
                           "The failed checks are listed below. Please send a new sample or re-sequence.", BIG))
        s.append(grid([["QC check", "Value", "Result"]] + gate_rows, [70, 60, 50]))
    else:
        s.append(Paragraph("5. Interpretation", H2))
        s.append(kv_table([("TB drug-resistance category", v11.get("category", "—")),
                           ("Mixed infection", f"{v10.get('verdict', '—')}" +
                            (f" (estimated second strain {v10.get('share')})" if v10.get("mixed") else "")),
                           ("Other pathogens", v02.get("others", "none") if v02.get("others") != "none" else "None detected"),
                           ("Points to note", v11.get("caveats", "none"))], cols=1))
        s.append(Paragraph("6. Drug susceptibility profile (genotypic)", H2))
        body = [["Drug", "Resistance mutation", "Mutation(s) and frequency", "WHO confidence"]]
        grp = None
        for r in results:
            if r["group"] != grp:
                grp = r["group"]
                body.append([f"<b>{esc(grp)}</b>", "", "", ""])
            muts = r["detail"] if r["call"] != "NOT DETECTED" else "—"
            body.append([f"{r['drug']} ({r['code']})", r["call"], muts, r["confidence"] or "—"])
        t = grid(body, [38, 30, 82, 30], badge_col=1)
        hdr_rows = [i for i, row in enumerate(body) if i and row[1] == ""]
        t.setStyle(TableStyle([("SPAN", (0, i), (-1, i)) for i in hdr_rows] +
                              [("BACKGROUND", (0, i), (-1, i), colors.HexColor("#E8EDF3")) for i in hdr_rows]))
        s.append(t)
        antibiogram(results, "antibiogram.png")
        iw, ih = ImageReader("antibiogram.png").getSize()
        s.append(Paragraph("7. Antibiogram", H2))
        s.append(Image("antibiogram.png", width=82 * mm, height=82 * mm * ih / iw))

        s.append(PageBreak())
        s.append(Paragraph("8. Variants used to call resistance", H2))
        rows = [["Position", "Ref", "Alt", "Type", "Freq %", "Depth", "Change", "Locus", "Gene", "Drug", "WHO grade"]]
        grouped = {}
        for e in evidence:
            grouped.setdefault((e["pos"], e["variant"]), []).append(e)
        for (pos, name), es in sorted(grouped.items(), key=lambda x: int(x[0][0])):
            v = variants.get(pos, {})
            gene = v.get("gene") or name.split("_")[0].split(" ")[0]
            ref, alt = v.get("ref", ""), v.get("alt", "")
            typ = "SNP" if len(ref) == len(alt) == 1 else ("MNP" if len(ref) == len(alt) else
                                                           ("INS" if len(alt) > len(ref) else "DEL")) if v else "DEL"
            change = name.split("_", 1)[1] if "_" in name else name
            rows.append([pos, ref or "—", alt or "—", typ, f"{100 * float(es[0]['af']):.1f}", v.get("dp", "—"),
                         change, locus.get(gene, "—"), gene, ", ".join(sorted({e["drug"] for e in es})),
                         "; ".join(sorted({e["grade"].split(") ")[-1] for e in es}))])
        s.append(grid(rows, [16, 9, 9, 11, 12, 13, 29, 15, 14, 22, 30]) if len(rows) > 1 else
                 Paragraph("No resistance-associated variants detected.", B))
        unc = [r for r in results if r["uncertain"]]
        if unc:
            s.append(Paragraph("Variants of uncertain significance (do not change the result)", H2))
            s.append(grid([["Drug", "Variants"]] + [[f"{r['drug']} ({r['code']})", r["uncertain"]] for r in unc],
                          [40, 140]))

    s.append(Paragraph("9. Definitions and criteria", H2))
    for txt in (
        "<b>Resistance mutation Detected</b>: a mutation graded by the WHO catalogue (2nd ed., 2023) as "
        "'Associated with resistance' (grade 1) or 'Associated with resistance - interim' (grade 2) was found "
        "at >= 10% of reads. <b>Not detected</b>: no such mutation found and all key positions were read at "
        ">= 10x. <b>Indeterminate</b>: no resistance mutation found, but key positions were not read well enough "
        "to exclude one.",
        "<b>WHO confidence</b>: grade assigned in the WHO catalogue. PAS is not graded by WHO; its result uses a "
        "literature list and is lower evidence. Pretomanid is inferred only from the WHO delamanid-pretomanid "
        "cross-resistance rule.",
        "<b>Categories (WHO 2021)</b>: Hr-TB: isoniazid-resistant, rifampicin-susceptible. RR-TB: rifampicin-"
        "resistant. MDR-TB: rifampicin and isoniazid. Pre-XDR-TB: RR/MDR-TB plus any fluoroquinolone. XDR-TB: "
        "pre-XDR-TB plus bedaquiline and/or linezolid.",
        "<b>Mixed infection</b>: >= 10 genome positions outside repetitive regions where a second allele is present "
        "at 10-90% with >= 5 reads, or more than one lineage detected. <b>Other pathogens</b>: organisms other than "
        "the MTB complex detected by Sylph at >= 95% genome identity and >= 1% abundance.",
    ):
        s.append(Paragraph(txt, B))
    s.append(Paragraph("10. Limitations and disclaimer", H2))
    s.append(Paragraph(
        "This is a genotypic prediction from whole-genome sequencing, not a phenotypic drug-susceptibility test. "
        "Resistance caused by mutations not yet graded by WHO, by mechanisms outside the analysed genes, or present "
        "in < 10% of the bacteria may not be detected. Where key positions were read 10-50 times, low-level "
        "(minority) resistance cannot be excluded. Results must be interpreted with clinical findings; treatment "
        "decisions follow national guidelines and the treating clinician's judgement.", B))
    s.append(Spacer(1, 10))
    s.append(kv_table([("Analysed by", a.operator), ("Reviewed by", ""), ("Authorised by", ""), ("Date", "")]))

    NumberedCanvas.footer = f"{a.lab_name}  |  {sid}  |  MTB_PIPELINE {a.version}  |  WHO catalogue 2023"
    out = os.path.join(a.out_dir, f"{sid}.clinical_report.pdf")
    SimpleDocTemplate(out, pagesize=A4, leftMargin=15 * mm, rightMargin=15 * mm, topMargin=12 * mm,
                      bottomMargin=15 * mm, title=f"{sid} TB WGS report", author=a.lab_name).build(
        s, canvasmaker=NumberedCanvas)

    w = PdfWriter()
    for p in [out] + sorted(a.step_pdfs, key=os.path.basename):
        w.append(p)
    w.write(os.path.join(a.out_dir, f"{sid}.all_steps.pdf"))

    json.dump({"schema": "mtb-step-summary/1", "sample_id": sid, "run_id": run, "step": "12_report",
               "status": "FAIL" if qc_fail else ("WARN" if qc_warn else "PASS"), "pipeline_version": a.version,
               "generated": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
               "values": {"qc": qc, "category": v11.get("category"), "mtb": mtb,
                          "detected": [r["code"] for r in results if r["call"] == "DETECTED"],
                          "indeterminate": [r["code"] for r in results if r["call"] == "INDETERMINATE"]},
               "metrics": [], "messages": [], "tools": [], "files": {}},
              open(os.path.join(a.out_dir, "12_report.summary.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
