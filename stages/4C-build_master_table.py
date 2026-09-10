#!/usr/bin/env python3
"""
Master table. One CSV, every species, every cholinesterase-family ortholog.

One row per SEQUENCE, not per species. Species with an endpoint and species
without are both present, distinguished by `set_membership`, so the file is the
single source of truth for the whole study and nothing downstream has to join
five tables again.

Non-AChE members (BChE, carboxylesterase, non-catalytic) are retained rather
than filtered. They are what the AChE calls were discriminated against, a
reviewer will want to see them, and the carboxylesterase counts per species are
themselves relevant to organophosphate sensitivity through sequestration.

  module purge
  python3 4C-build_master_table.py --classify CLASSIFY --tree TREE --odb ODB \\
      --final FINAL --final-set FINAL_SET --ecotox ../supp/Chlorpyrifos_cleaned_gmean.txt \\
      --out master_table.csv
"""

import argparse, csv, math, os, sys
from collections import defaultdict

CLADE_LABEL = {"AChE_vert": "Vertebrate AChE", "AChE_insect": "Arthropod AChE",
               "AChE_inv": "Invertebrate AChE", "BChE": "BChE",
               "CCE": "Carboxylesterase", "dead": "Non-catalytic"}
CLASS_ALIAS = {"Actinopteri": "Actinopterygii", "Maxillopoda": "Copepoda",
               "Hexanauplia": "Copepoda", "Turbellaria": "Rhabditophora"}
RESIDUES = ["Tyr70", "Asp72", "Trp84", "Gly118", "Gly119", "Tyr121", "Glu199",
            "Ser200", "Ala201", "Trp279", "Phe288", "Phe290", "Glu327",
            "Phe330", "Tyr334", "His440"]


def load(path, required=True):
    if not os.path.exists(path):
        if required:
            sys.exit(f"missing: {path}")
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classify", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--odb", required=True)
    ap.add_argument("--final", required=True)
    ap.add_argument("--final-set", required=True)
    ap.add_argument("--ecotox", required=True)
    ap.add_argument("--out", default="master_table.csv")
    a = ap.parse_args()

    calls = load(os.path.join(a.classify, "active_site_calls.tsv"))
    tassign = {r["sequence_id"]: r
               for r in load(os.path.join(a.tree, "tree_assignments.tsv"),
                             required=False)}
    reps = {r["sequence_id"]: r
            for r in load(os.path.join(a.tree, "representatives.tsv"))}
    inv = load(os.path.join(a.odb, "species_inventory.tsv"))
    cov = {r["Species_Name"]: r
           for r in load(os.path.join(a.odb, "ecotox_coverage.tsv"),
                         required=False)}
    final = load(os.path.join(a.final_set, "endpoint_final.tsv"))
    status = {r["Species_Name"]: r
              for r in load(os.path.join(a.final, "endpoint_species.tsv"))}
    eco = load(a.ecotox)

    tax = {}
    for r in inv:
        k = (r.get("binomial") or r.get("organism_name", "")).replace(" ", "_")
        tax.setdefault(k, r)

    # endpoint, keyed by the ODB species tag that supplied the sequence
    ep_by_tag = defaultdict(list)
    for r in final:
        ep_by_tag[r["odb_species"]].append(r)
    eco_by_name = {r["Species_Name"].strip(): r for r in eco}

    rows = []
    for c in calls:
        sid, sp = c["sequence_id"], c["species"]
        t = tassign.get(sid, {})
        tx = tax.get(sp, {})
        grp = t.get("anchor_group", "")
        eps = ep_by_tag.get(sp, [])
        ep = eps[0] if eps else None
        raw = eco_by_name.get(ep["Species_Name"], {}) if ep else {}

        # OrthoDB header: species__geneid__orthodb__ogid
        parts = sid.split("__")
        gene_id = parts[1] if len(parts) > 1 else ""
        og_id = parts[3] if len(parts) > 3 else ""

        row = {
            "sequence_id": sid,
            "odb_species": sp.replace("_", " "),
            "odb_gene_id": gene_id,
            "orthogroup": og_id,
            "ncbi_taxid": tx.get("ncbi_taxid", ""),
            "class": CLASS_ALIAS.get(tx.get("class", ""), tx.get("class", "")),
            "order": tx.get("order", ""),
            "family": tx.get("family", ""),
            # --- study membership ---
            "set_membership": ("endpoint" if eps else "prediction"),
            "ecotox_species": "; ".join(x["Species_Name"] for x in eps),
            "log10_LC50": ep["log10_LC50"] if ep else "",
            "LC50_mg_L": (f"{10**float(ep['log10_LC50']):.4g}" if ep else ""),
            "conc_type": raw.get("Conc_Type", ""),
            "pub_year": raw.get("Pub_Year", ""),
            "sequence_source": ("direct" if ep and ep.get("odb_match") == "exact"
                                else ("congener proxy" if ep else "")),
            "tier": ep.get("tier", "") if ep else "",
            # --- classification ---
            "residue_call": c["call"],
            "residue_confidence": c["confidence"],
            "tree_clade": CLADE_LABEL.get(grp, grp),
            "tree_clade_raw": grp,
            "clade_anchors": t.get("clade_anchors", ""),
            "clade_unambiguous": t.get("clade_unambiguous", ""),
            "tree_vs_residue": t.get("tree_vs_residue", ""),
            "paralogue": c.get("paralogue", ""),
            "pid_ace1": c.get("pid_ace1", ""),
            "pid_ace2": c.get("pid_ace2", ""),
            "is_representative": sid in reps,
            "selection_basis": reps.get(sid, {}).get("selection_basis", ""),
            # --- active site ---
            "triad_intact": c.get("triad_intact", ""),
            "n_gorge_aromatic": c.get("n_gorge_aromatic", ""),
            "n_gaps_at_sites": c.get("n_gaps_at_sites", ""),
            "active_site_string": "".join(c.get(p, "-") for p in RESIDUES),
        }
        for p in RESIDUES:
            row[f"res_{p}"] = c.get(p, "")
        rows.append(row)

    # ECOTOX species that never reached the classifier at all
    seen_eco = {x for r in rows for x in r["ecotox_species"].split("; ") if x}
    for name, r in eco_by_name.items():
        if name in seen_eco:
            continue
        s = status.get(name, {})
        try:
            y = round(math.log10(float(r["Conc_Mean_Std"])), 4)
        except (ValueError, KeyError, TypeError):
            y = ""
        blank = {k: "" for k in rows[0]} if rows else {}
        blank.update({
            "sequence_id": "", "odb_species": "", "class": r.get("Class", ""),
            "set_membership": "endpoint_no_sequence",
            "ecotox_species": name, "log10_LC50": y,
            "LC50_mg_L": r.get("Conc_Mean_Std", ""),
            "conc_type": r.get("Conc_Type", ""), "pub_year": r.get("Pub_Year", ""),
            "sequence_source": "none",
            "selection_basis": s.get("status", "no_sequence"),
            "is_representative": False,
        })
        rows.append(blank)

    cols = list(rows[0])
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    print(f"{len(rows)} rows -> {a.out}")
    print("\nby set_membership:")
    for k, v in Counter(r["set_membership"] for r in rows).most_common():
        print(f"  {k}: {v}")
    print("\nby tree clade:")
    for k, v in Counter(r["tree_clade"] for r in rows if r["tree_clade"]).most_common():
        print(f"  {k}: {v}")
    print(f"\nrepresentatives: {sum(1 for r in rows if r['is_representative'])}")
    print(f"species covered:  "
          f"{len({r['odb_species'] for r in rows if r['odb_species']})}")
    print("\nOne row per sequence. Filter is_representative == True for the "
          "modelling set; the rest are the family members those calls were "
          "discriminated against, and the carboxylesterase counts per species "
          "are worth keeping since sequestration is a competing explanation "
          "for interspecies OP sensitivity.")


if __name__ == "__main__":
    main()
