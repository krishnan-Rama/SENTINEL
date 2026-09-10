#!/usr/bin/env python3
"""
Build the final analysis set.

Joins the tree representatives back to the ECOTOX endpoints and the
supplementary species, and emits the two things everything downstream needs:

  endpoint_final.tsv / .faa    species with an LC50 and one chosen AChE.
                               This is your n.
  prediction_final.tsv / .faa  species with an AChE and no LC50. The SSD
                               prediction set.

Confidence tiers, carried as a covariate rather than used as a filter:
  A  direct species match in OrthoDB, AChE in the expected clade, high
     residue confidence, single candidate
  B  as A but a genus proxy, or medium residue confidence, or several
     candidates in the clade
  C  anything that needed a judgement call. Fit the model with and without.

  module purge
  python3 4B-build_final_set.py --final FINAL --tree TREE --outdir FINAL_SET
"""

import argparse, csv, os, sys
from collections import Counter, defaultdict


def load(path):
    if not os.path.exists(path):
        sys.exit(f"missing: {path}")
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path):
    recs, h, s = [], None, []
    for l in open(path):
        if l.startswith(">"):
            if h: recs.append((h, "".join(s)))
            h, s = l[1:].rstrip("\n"), []
        elif l.strip():
            s.append(l.strip())
    if h: recs.append((h, "".join(s)))
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final", required=True, help="FINAL/ from 4A-finalise_endpoints.py")
    ap.add_argument("--tree", required=True, help="TREE/ from 3B-assign_by_tree.py")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    eps = load(os.path.join(a.final, "endpoint_species.tsv"))
    reps = load(os.path.join(a.tree, "representatives.tsv"))
    seqs = {h.split()[0]: s
            for h, s in read_fasta(os.path.join(a.tree, "representatives.faa"))}

    rep_by_sp = {r["species"]: r for r in reps}

    def tier(ep, rp):
        if rp.get("selection_basis", "").upper().startswith("NO ACHE"):
            return "C", "no AChE clade placement"
        why = []
        t = "A"
        if ep.get("odb_match") != "exact":
            t, w = "B", "genus proxy, not direct species evidence"
            why.append(w)
        if rp.get("residue_conf") != "high":
            t = "B" if t == "A" else t
            why.append(f"residue confidence {rp.get('residue_conf')}")
        if rp.get("tree_vs_residue") == "DISAGREE":
            t = "B" if t == "A" else t
            why.append("tree and residue call differ")
        if str(rp.get("clade_unambiguous", "")).lower() == "false":
            t = "C"
            why.append("clade contains more than one anchor group")
        try:
            if int(rp.get("n_candidates", "1")) > 1:
                t = "B" if t == "A" else t
                why.append(f"{rp.get('n_candidates')} candidates in clade")
        except ValueError:
            pass
        return t, "; ".join(why) or "direct, unambiguous"

    rows, missing = [], []
    for ep in eps:
        tags = [t for t in (ep.get("odb_tags") or "").split(";") if t]
        rp = next((rep_by_sp[t] for t in tags if t in rep_by_sp), None)
        if rp is None:
            missing.append((ep["Species_Name"], ep.get("ecotox_class", ""),
                            ep.get("status", "")))
            continue
        t, why = tier(ep, rp)
        rows.append({
            "Species_Name": ep["Species_Name"],
            "ecotox_class": ep.get("ecotox_class", ""),
            "log10_LC50": ep["log10_LC50"],
            "odb_match": ep.get("odb_match", ""),
            "odb_species": rp["species"],
            "sequence_id": rp["sequence_id"],
            "clade_group": rp.get("anchor_group", ""),
            "clade_anchors": rp.get("clade_anchors", ""),
            "residue_call": rp.get("residue_call", ""),
            "residue_conf": rp.get("residue_conf", ""),
            "paralogue": rp.get("paralogue", ""),
            "suggested_paralogue": ep.get("suggested_paralogue", ""),
            "n_candidates": rp.get("n_candidates", ""),
            "tier": t, "tier_reason": why,
        })

    if not rows:
        sys.exit("no endpoint species could be joined to a representative")

    with open(os.path.join(a.outdir, "endpoint_final.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    n_seq = 0
    with open(os.path.join(a.outdir, "endpoint_final.faa"), "w") as fh:
        for r in rows:
            s = seqs.get(r["sequence_id"])
            if not s:
                continue
            fh.write(f">{r['Species_Name'].replace(' ', '_')}|{r['tier']}"
                     f"|{r['clade_group']}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i+60] + "\n")
            n_seq += 1

    ys = [float(r["log10_LC50"]) for r in rows]
    print(f"== ENDPOINT SET: n = {len(rows)} ==")
    print(f"  log10 LC50 {min(ys):.2f} to {max(ys):.2f} "
          f"(span {max(ys)-min(ys):.2f}; full dataset 5.76)")
    print(f"  sequences written: {n_seq}")
    print("\n  by class:")
    for c, n in Counter(r["ecotox_class"] for r in rows).most_common():
        sub = [float(x["log10_LC50"]) for x in rows if x["ecotox_class"] == c]
        print(f"    {c:16s} {n:3d}   {min(sub):6.2f} to {max(sub):6.2f}")
    print("\n  by tier:")
    for t, n in sorted(Counter(r["tier"] for r in rows).items()):
        print(f"    {t}: {n}")
    ta = [r for r in rows if r["tier"] == "A"]
    if ta:
        ya = [float(r["log10_LC50"]) for r in ta]
        print(f"  tier A alone: n={len(ta)}, span {max(ya)-min(ya):.2f}")

    if missing:
        print(f"\n  endpoint species with no representative: {len(missing)}")
        for sp, cls, st in missing[:12]:
            print(f"    {sp:30s} {cls:16s} {st}")
        if len(missing) > 12:
            print(f"    ... and {len(missing)-12} more")

    # ---- prediction set ----------------------------------------------------
    ep_tags = {r["odb_species"] for r in rows}
    pred = [r for r in reps if r["species"] not in ep_tags
            and not r.get("selection_basis", "").upper().startswith("NO ACHE")]
    with open(os.path.join(a.outdir, "prediction_final.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(pred[0]))
        w.writeheader(); w.writerows(pred)
    with open(os.path.join(a.outdir, "prediction_final.faa"), "w") as fh:
        for r in pred:
            s = seqs.get(r["sequence_id"])
            if not s:
                continue
            fh.write(f">{r['species']}|pred|{r.get('anchor_group','')}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i+60] + "\n")
    print(f"\n== PREDICTION SET: {len(pred)} species ==")
    for c, n in Counter(r.get("class", "") for r in pred).most_common(12):
        print(f"    {c or 'unknown':16s} {n}")

    print(f"\nwritten to {a.outdir}/")
    print("\nThe endpoint set is what the model is fitted and tested on. The "
          "prediction set has no LC50, so nothing about it can be validated "
          "within this study. Report the two separately and never pool them "
          "into a single performance figure.")


if __name__ == "__main__":
    main()
