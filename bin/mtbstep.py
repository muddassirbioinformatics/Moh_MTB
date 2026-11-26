"""Shared helpers for step summary scripts (schema mtb-step-summary/1, see step_report.py)."""
import datetime, json
import yaml

# chart tokens shared by all step plots (reference palette, light print surface)
SERIES1, SERIES2 = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"


class Step:
    def __init__(self, step, sample_id, run_id, version, thresholds_file):
        th = yaml.safe_load(open(thresholds_file))
        self.th_version = th.get("version")
        self.t = th.get(step, {})
        self.prov = self.t.get("provisional", False)
        self.d = {"schema": "mtb-step-summary/1", "sample_id": sample_id, "run_id": run_id, "step": step,
                  "pipeline_version": version, "values": {}, "metrics": [], "messages": [], "tools": [],
                  "files": {"inputs": [], "outputs": []}}

    def metric(self, name, plain, value, unit="", threshold="reported only", status="INFO"):
        self.d["metrics"].append({"name": name, "plain": plain, "value": value, "unit": unit,
                                  "threshold": threshold, "status": status,
                                  "provisional": self.prov and status != "INFO"})

    def msg(self, text):
        self.d["messages"].append(text)

    def tool(self, name, version, container, command=None):
        self.d["tools"].append({"name": name, "version": version, "container": container, "command": command})

    def output(self, path):
        self.d["files"]["outputs"].append({"path": path})

    def write(self, out, **values):
        st = [m["status"] for m in self.d["metrics"]]
        self.d["status"] = "FAIL" if "FAIL" in st else "WARN" if "WARN" in st else "PASS"
        self.d["values"].update(values)
        self.d["generated"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        self.d["thresholds_version"] = self.th_version
        json.dump(self.d, open(out, "w"), indent=2)


def grade(value, fail_below=None, warn_below=None):
    if fail_below is not None and value < fail_below:
        return "FAIL"
    if warn_below is not None and value < warn_below:
        return "WARN"
    return "PASS"


def style_axes(ax):
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.tick_params(colors=INK2, labelsize=7.5, length=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)


def key_depth(depth, start0, end):
    """Depth used for a catalogue position: minimum for SNVs/short changes (<= 50 bp); median for long
    catalogue deletions, whose span can run past the gene (one low base there must not void the drug)."""
    seg = depth[start0:end]
    return int(seg.min()) if end - start0 <= 50 else int(__import__("numpy").median(seg))


def read_text(path):
    return open(path).read().strip()
