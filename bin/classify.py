#!/usr/bin/env python3
"""Step 11: WHO programmatic resistance category from the step 09 drug results.

WHO definitions (2021): Hr-TB = INH-R, RIF-S; RR-TB = RIF-R; MDR-TB = RIF-R + INH-R;
pre-XDR = RR/MDR + resistance to any fluoroquinolone (LFX or MFX);
XDR = pre-XDR + resistance to BDQ and/or LZD.
An Indeterminate drug that could change the category is stated, never assumed susceptible.
"""
import argparse, csv, json
from mtbstep import Step


def main():
    ap = argparse.ArgumentParser()
    for k in ("drug_results", "mixed", "sample_id", "run_id", "thresholds", "version", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    s = Step("11_classification", a.sample_id, a.run_id, a.version, a.thresholds)
    res = {r["code"]: r for r in csv.DictReader(open(a.drug_results), delimiter="\t")}
    R = lambda c: res[c]["call"] == "DETECTED"
    I = lambda c: res[c]["call"] == "INDETERMINATE"
    fq_r, fq_i = R("LFX") or R("MFX"), (I("LFX") or I("MFX")) and not (R("LFX") or R("MFX"))
    grpA_r, grpA_i = R("BDQ") or R("LZD"), (I("BDQ") or I("LZD")) and not (R("BDQ") or R("LZD"))
    caveats = []

    if R("RIF"):
        cat = "MDR-TB" if R("INH") else "RR-TB"
        if I("INH"):
            caveats.append("isoniazid Indeterminate: MDR-TB cannot be excluded")
        if fq_r:
            cat = "Pre-XDR-TB"
            if grpA_r:
                cat = "XDR-TB"
            elif grpA_i:
                caveats.append("bedaquiline/linezolid Indeterminate: XDR-TB cannot be excluded")
        elif fq_i:
            caveats.append("fluoroquinolone Indeterminate: pre-XDR-TB cannot be excluded")
    elif I("RIF"):
        cat = "Hr-TB (rifampicin Indeterminate)" if R("INH") else "Indeterminate"
        caveats.append("rifampicin Indeterminate: RR-TB cannot be excluded")
    elif R("INH"):
        cat = "Hr-TB"
    else:
        other = [res[c]["drug"] for c in res if R(c)]
        cat = ("Other drug resistance (" + ", ".join(other) + ")") if other else \
              "No resistance detected (to the drugs tested)"
        if I("INH"):
            caveats.append("isoniazid Indeterminate: Hr-TB cannot be excluded")

    minority = [res[c]["drug"] for c in res if R(c) and "minority population" in res[c]["detail"]]
    if minority:
        caveats.append("based partly on minority (heteroresistant) populations for: " + ", ".join(minority))
    mixed = json.load(open(a.mixed))["values"].get("mixed")
    if mixed:
        caveats.append("mixed-strain infection detected (step 10)")
    resistant = [f"{res[c]['drug']} ({c})" for c in res if R(c)]

    s.metric("WHO resistance category", "Type of drug-resistant TB", cat, threshold="WHO 2021 definitions",
             status="WARN" if caveats else "PASS")
    s.metric("Resistant to", "Drugs with resistance detected", ", ".join(resistant) or "none")
    s.metric("Caveats", "Points to note", "; ".join(caveats) or "none")
    s.msg("Definitions (WHO 2021): Hr-TB = isoniazid-resistant, rifampicin-susceptible; RR-TB = rifampicin-"
          "resistant; MDR-TB = rifampicin + isoniazid; pre-XDR-TB = RR/MDR-TB + any fluoroquinolone; XDR-TB = "
          "pre-XDR-TB + bedaquiline and/or linezolid.")
    s.msg("This is a genotypic classification. Treatment decisions follow national guidelines and the treating "
          "clinician's judgement.")
    s.tool("classify.py", a.version, "mtb-report:1.1.0")
    s.output("11_classification.summary.json")
    s.write(a.out, category=cat, caveats="; ".join(caveats) or "none", resistant=", ".join(resistant) or "none")


if __name__ == "__main__":
    main()
