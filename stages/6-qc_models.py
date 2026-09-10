#!/usr/bin/env python3
"""
QC the ColabFold models where it matters: the active site.

Whole-model pLDDT is the wrong statistic for this study. A model can average 92
while the acyl pocket is ragged, and it will still dock and still return a
number. What decides whether a species is usable is the confidence at the 16
Torpedo-numbered positions your classifier reads and your docking box sits on.

Two things are measured:

  1. Per-residue pLDDT at the active-site positions, mapped through the same
     alignment frame used for trimming, so the same structural position is
     inspected in every species.

  2. Between-seed active-site RMSD. All five models are superposed on rank 1 by
     their CA atoms, then the RMSD of the active-site CAs is computed. This is
     your NOISE FLOOR. Any between-species difference in gorge geometry smaller
     than this number is not measurable with these structures, and knowing it
     before docking is the difference between a result and an artefact.

  module purge
  module load Python/3.11.5-GCCcore-13.2.0 SciPy-bundle/2025.07-gfbf-2025b
  python3 6-qc_models.py --models FOLD/models --aln ../CLASSIFY/aligned.a2m \\
      --final-set ../FINAL_SET --outdir QC
"""

import argparse, csv, glob, json, math, os, re, shutil, sys
from collections import defaultdict

try:
    import numpy as np
except ImportError:
    sys.exit("numpy required. module load SciPy-bundle/2025.07-gfbf-2025b")

RESIDUES = ["Tyr70", "Asp72", "Trp84", "Gly118", "Gly119", "Tyr121", "Glu199",
            "Ser200", "Ala201", "Trp279", "Phe288", "Phe290", "Glu327",
            "Phe330", "Tyr334", "His440"]
ELBOW = re.compile(r"G[EQ]S[AG]G")
# The acyl pocket and the catalytic serine are what the docking box and the
# classifier depend on. A model can be excellent elsewhere and useless here.
CRITICAL = ["Ser200", "Glu327", "His440", "Phe288", "Phe290", "Trp84"]


def read_fasta(p):
    recs, h, s = [], None, []
    for l in open(p):
        if l.startswith(">"):
            if h: recs.append((h, "".join(s)))
            h, s = l[1:].rstrip("\n"), []
        elif l.strip(): s.append(l.strip())
    if h: recs.append((h, "".join(s)))
    return recs


def kabsch_transform(P, Q):
    """Rotation and centroids that superpose P onto Q using all given atoms."""
    pc, qc = P.mean(0), Q.mean(0)
    V, S, Wt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return R, pc, qc


def kabsch_rmsd(P, Q, idx=None):
    """Superpose P onto Q using ALL atoms, then RMSD over idx only."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    Pr = Pc @ R
    sel = slice(None) if idx is None else idx
    diff = Pr[sel] - Qc[sel]
    return float(np.sqrt((diff ** 2).sum() / len(diff)))


def parse_pdb_atoms(path, resnums):
    """(resnum, atomname) -> coord, for the given residues, heavy atoms only."""
    want, out = set(resnums), {}
    for l in open(path):
        if not l.startswith("ATOM"):
            continue
        r = int(l[22:26])
        if r not in want:
            continue
        nm = l[12:16].strip()
        if nm.startswith("H") or l[76:78].strip() == "H":
            continue
        out[(r, nm)] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return out


def parse_pdb_ca(path):
    """residue number -> CA coordinate, plus resname."""
    ca, name = {}, {}
    for l in open(path):
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            r = int(l[22:26])
            ca[r] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            name[r] = l[17:20].strip()
    return ca, name


def parse_pdb_atom(path, resnum, atom):
    for l in open(path):
        if l.startswith("ATOM") and int(l[22:26]) == resnum \
                and l[12:16].strip() == atom:
            return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", required=True)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--final-set", default=None)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--end", type=int, default=537)
    ap.add_argument("--min-plddt", type=float, default=70.0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    # ---- rebuild the trimming frame, exactly as prepare_for_folding did ----
    aln = {}
    for h, s in read_fasta(a.aln):
        aln.setdefault(h.split()[0], s)
    tor = next((k for k in aln if "Torpedo" in k), None)
    if tor is None:
        sys.exit("no Torpedo anchor in the alignment")

    res2col, r, col = {}, 0, 0
    for ch in aln[tor]:
        if ch == "-": col += 1
        elif ch == ".": continue
        elif ch.isupper(): r += 1; res2col[r] = col; col += 1
        else: r += 1; res2col[r] = None
    raw = "".join(c for c in aln[tor] if c not in "-.").upper()
    m = ELBOW.search(raw)
    if not m:
        sys.exit("elbow motif not found")
    offset = (m.start() + 3) - 200

    def col_of(mature):
        up = mature + offset
        c = res2col.get(up)
        while c is None and up in res2col:
            up += 1; c = res2col.get(up)
        return c

    c_start, c_end = col_of(a.start), col_of(a.end)
    site_cols = {lab: col_of(int(re.search(r"\d+", lab).group()))
                 for lab in RESIDUES}
    print(f"offset {offset}; trimmed columns {c_start}-{c_end}")

    def site_positions(seq_id):
        """active-site label -> 1-based index in the TRIMMED sequence."""
        if seq_id not in aln:
            return None
        out, kept, cur = {}, 0, 0
        want = {v: k for k, v in site_cols.items()}
        for ch in aln[seq_id]:
            if ch == ".":
                continue
            if ch == "-":
                cur += 1
                continue
            inside = c_start <= cur <= c_end
            if inside:
                kept += 1
            if ch.isupper() and cur in want and inside:
                out[want[cur]] = kept
            if ch.isupper():
                cur += 1
        return out

    # ---- iterate models ---------------------------------------------------
    best_dir = os.path.join(a.outdir, "best_models")
    os.makedirs(best_dir, exist_ok=True)
    rows = []
    pdbs = sorted(glob.glob(os.path.join(a.models, "**",
                                         "*_relaxed_rank_001_*.pdb"),
                            recursive=True))
    print(f"{len(pdbs)} rank-1 models found")

    for p1 in pdbs:
        base = os.path.basename(p1).split("_relaxed_rank_001_")[0]
        d = os.path.dirname(p1)
        sc = sorted(glob.glob(os.path.join(d, base + "_scores_rank_001_*.json")))
        if not sc:
            continue
        j = json.load(open(sc[0]))
        plddt = np.array(j.get("plddt", []), dtype=float)
        ptm = j.get("ptm", "")

        sp = base.split("__")[0]
        pos = site_positions(base)
        row = {"species": sp, "sequence_id": base, "n_residues": len(plddt),
               "ptm": round(float(ptm), 4) if ptm != "" else "",
               "plddt_mean": round(float(plddt.mean()), 2) if len(plddt) else ""}

        if not pos:
            row.update({"verdict": "no_alignment_mapping"})
            rows.append(row); continue

        vals, gapped = {}, []
        for lab in RESIDUES:
            i = pos.get(lab)
            if i is None or i > len(plddt):
                gapped.append(lab); continue
            vals[lab] = float(plddt[i - 1])
        row["n_sites_gapped"] = len(gapped)
        row["sites_gapped"] = ";".join(gapped)
        if vals:
            row["plddt_site_mean"] = round(sum(vals.values()) / len(vals), 2)
            row["plddt_site_min"] = round(min(vals.values()), 2)
        crit = [vals[c] for c in CRITICAL if c in vals]
        row["plddt_critical_min"] = round(min(crit), 2) if crit else ""
        for lab in RESIDUES:
            row[f"plddt_{lab}"] = round(vals[lab], 1) if lab in vals else ""

        # ---- between-seed active-site RMSD --------------------------------
        seeds = sorted(glob.glob(os.path.join(d, base + "_relaxed_rank_00*_*.pdb")))
        ref_ca, _ = parse_pdb_ca(p1)
        common = sorted(ref_ca)
        site_idx = [common.index(pos[l]) for l in RESIDUES
                    if l in pos and pos[l] in ref_ca]
        rms = []
        Q = np.array([ref_ca[k] for k in common])
        for s in seeds:
            if s == p1:
                continue
            ca, _ = parse_pdb_ca(s)
            if set(ca) != set(ref_ca):
                continue
            P = np.array([ca[k] for k in common])
            if site_idx:
                rms.append(kabsch_rmsd(P, Q, idx=site_idx))
        row["n_seeds_compared"] = len(rms)
        row["seed_rmsd_site_mean"] = round(float(np.mean(rms)), 3) if rms else ""
        row["seed_rmsd_site_max"] = round(float(np.max(rms)), 3) if rms else ""

        # ALL-ATOM active-site RMSD. Backbone is not what determines gorge
        # volume, acyl-pocket shape or near-attack distance; side-chain
        # rotamers are, and they scatter between seeds far more than CA does.
        # A CA-only noise floor therefore understates the uncertainty on
        # exactly the quantities this study claims species differences in.
        site_res = [pos[l] for l in RESIDUES if l in pos and pos[l] in ref_ca]
        aa = []
        if site_res:
            ref_at = parse_pdb_atoms(p1, site_res)
            for sfile in seeds:
                if sfile == p1:
                    continue
                ca, _ = parse_pdb_ca(sfile)
                if set(ca) != set(ref_ca):
                    continue
                at = parse_pdb_atoms(sfile, site_res)
                keys = sorted(set(ref_at) & set(at))
                if not keys:
                    continue
                P = np.array([ca[k] for k in common])
                R, pc, qc = kabsch_transform(P, Q)   # superpose on ALL CA
                A = (np.array([at[k] for k in keys]) - pc) @ R
                B = np.array([ref_at[k] for k in keys]) - qc
                aa.append(float(np.sqrt(((A - B) ** 2).sum() / len(keys))))
        row["n_site_atoms"] = len(set(ref_at)) if site_res else 0
        row["seed_rmsd_site_allatom_mean"] = (round(float(np.mean(aa)), 3)
                                              if aa else "")
        row["seed_rmsd_site_allatom_max"] = (round(float(np.max(aa)), 3)
                                             if aa else "")

        # catalytic Ser Og: the docking box centre
        if "Ser200" in pos:
            og = parse_pdb_atom(p1, pos["Ser200"], "OG")
            if og:
                row["ser_resnum"] = pos["Ser200"]
                row["ser_og_x"], row["ser_og_y"], row["ser_og_z"] = \
                    [round(v, 3) for v in og]

        pm = row.get("plddt_critical_min", "")
        if gapped and any(g in CRITICAL for g in gapped):
            row["verdict"] = "FAIL_critical_site_missing"
        elif pm == "" :
            row["verdict"] = "FAIL_no_site_plddt"
        elif pm < 50:
            row["verdict"] = "FAIL_very_low"
        elif pm < a.min_plddt:
            row["verdict"] = "REVIEW_low"
        else:
            row["verdict"] = "pass"
            shutil.copy(p1, os.path.join(best_dir, f"{sp}.pdb"))
        rows.append(row)

    cols = ["species", "sequence_id", "verdict", "n_residues", "ptm",
            "plddt_mean", "plddt_site_mean", "plddt_site_min",
            "plddt_critical_min", "n_sites_gapped", "sites_gapped",
            "n_seeds_compared", "seed_rmsd_site_mean", "seed_rmsd_site_max",
            "n_site_atoms", "seed_rmsd_site_allatom_mean",
            "seed_rmsd_site_allatom_max",
            "ser_resnum", "ser_og_x", "ser_og_y", "ser_og_z"] + \
           [f"plddt_{l}" for l in RESIDUES]
    with open(os.path.join(a.outdir, "model_qc.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols,
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    from collections import Counter
    print("\n== verdicts ==")
    for k, v in Counter(r["verdict"] for r in rows).most_common():
        print(f"  {k}: {v}")

    ok = [r for r in rows if r["verdict"] == "pass"]
    if ok:
        sm = [r["plddt_site_mean"] for r in ok if r.get("plddt_site_mean") != ""]
        cm = [r["plddt_critical_min"] for r in ok if r.get("plddt_critical_min") != ""]
        rr = [r["seed_rmsd_site_mean"] for r in ok
              if r.get("seed_rmsd_site_mean") not in ("", None)]
        print(f"\npassing models: {len(ok)}")
        print(f"  active-site pLDDT mean: {min(sm):.1f} to {max(sm):.1f} "
              f"(median {sorted(sm)[len(sm)//2]:.1f})")
        print(f"  critical-residue pLDDT min: {min(cm):.1f} to {max(cm):.1f}")
        aa = [r["seed_rmsd_site_allatom_mean"] for r in ok
              if r.get("seed_rmsd_site_allatom_mean") not in ("", None)]
        if rr:
            print(f"\n== between-seed active-site RMSD (the noise floor) ==")
            print(f"  backbone CA : {min(rr):.2f} to {max(rr):.2f} A, "
                  f"median {sorted(rr)[len(rr)//2]:.2f}")
        if aa:
            print(f"  ALL ATOM    : {min(aa):.2f} to {max(aa):.2f} A, "
                  f"median {sorted(aa)[len(aa)//2]:.2f}")
            print("\n  Use the ALL-ATOM figure. Gorge volume and near-attack "
                  "geometry are set by side-chain rotamers, not backbone, so "
                  "the CA number understates the uncertainty on every metric "
                  "you intend to compare across species. Any between-species "
                  "difference below the all-atom value is seed variance.")

    bad = [r for r in rows if r["verdict"] != "pass"]
    if bad:
        print(f"\nnot passing ({len(bad)}):")
        for r in sorted(bad, key=lambda x: str(x.get("plddt_critical_min"))):
            print(f"  {r['species']:32s} {r['verdict']:28s} "
                  f"crit_min={r.get('plddt_critical_min','')} "
                  f"gapped={r.get('sites_gapped','')}")

    print(f"\n  {a.outdir}/model_qc.tsv")
    print(f"  {best_dir}/  ({len(ok)} models, one per species)")
    print("\nser_og_* is the catalytic serine Og coordinate in each model. That "
          "is the docking box centre, and it is per-model: do not reuse one "
          "species' coordinates for another.")


if __name__ == "__main__":
    main()
