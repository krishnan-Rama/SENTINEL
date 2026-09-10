#!/usr/bin/env python3
"""
Post-docking analysis, in the order the conclusions depend on each other.

  1  AFFINITY NOISE. Three seeds per receptor-ligand pair were run for the same
     reason five models per species were: to find out whether the between-species
     spread in score exceeds the between-seed spread. ICC per ligand. If a
     ligand's ICC is low, its affinities are not a species measurement and no
     regression on them means anything.

  2  NEAR-ATTACK GEOMETRY, for the oxon only. A Vina score is a non-covalent
     affinity; organophosphate potency is k2/Kd and nothing in a docking score
     touches k2. What a pose CAN tell you is whether the reactive geometry is
     reachable: the distance from Ser Og to the phosphorus, and how close the
     Og-P-O(leaving) angle comes to the 180 degrees required for in-line
     phosphoryl transfer. These are physically interpretable and comparable
     across species in a way that a score is not.

  3  CONTACT FINGERPRINTS in Torpedo numbering. Every model residue is mapped
     through the alignment to a Torpedo mature position, so "the ligand touches
     position 288" means the same thing in a fish and a water flea. This is the
     interaction-level comparison, and it is the part that can show WHICH
     residues differ where sensitivity differs.

  4  CONTROL DIAGNOSTIC. Each ligand's affinity is correlated against log10
     LC50 separately. Chlorpyrifos itself cannot inhibit AChE without CYP
     desulfuration, and TCP cannot phosphorylate anything at all. If either
     predicts LC50 as well as the oxon does, the pipeline is ranking
     lipophilicity or generic gorge accommodation rather than target
     engagement. That is the single most important number in this script and
     it is designed to be able to falsify the whole approach.

  module load micromamba/2.8.1
  export PATH=/mnt/ecotox/GROUP-smbpk/c23048124/envs/dock/bin:$PATH
  python3 9A-analyse_docking.py --docked DOCKED --receptors RECEPTORS/prepped \\
      --prep RECEPTORS/receptor_prep.tsv --aln ../CLASSIFY/aligned.a2m \\
      --endpoints ../REPORT/analysis_set.tsv --outdir DOCKANALYSIS
"""

import argparse, csv, glob, math, os, re, sys
from collections import defaultdict, Counter

try:
    import numpy as np
except ImportError:
    sys.exit("numpy required")

ELBOW = re.compile(r"G[EQ]S[AG]G")
OXON = "chlorpyrifos_oxon"
CONTACT_CUT = 4.0
OXYANION = ["Gly118", "Gly119", "Ala201"]
SITES = ["Tyr70", "Asp72", "Trp84", "Gly118", "Gly119", "Tyr121", "Glu199",
         "Ser200", "Ala201", "Trp279", "Phe288", "Phe290", "Glu327",
         "Phe330", "Tyr334", "His440"]


def read_fasta(p):
    recs, h, s = [], None, []
    for l in open(p):
        if l.startswith(">"):
            if h: recs.append((h, "".join(s)))
            h, s = l[1:].rstrip("\n"), []
        elif l.strip(): s.append(l.strip())
    if h: recs.append((h, "".join(s)))
    return recs


def read_tsv(p):
    if not os.path.exists(p):
        return []
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pdb_heavy(path):
    xyz, resi, resn, name = [], [], [], []
    for l in open(path):
        if l.startswith("ATOM"):
            nm = l[12:16].strip()
            if nm.startswith("H"):
                continue
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            resi.append(int(l[22:26])); resn.append(l[17:20].strip()); name.append(nm)
    return np.array(xyz), np.array(resi), np.array(resn, dtype=object), \
        np.array(name, dtype=object)


def first_pose(pdbqt):
    """Best (first) MODEL: coords, autodock types, and the reported affinity."""
    xyz, typ, aff = [], [], None
    inmodel = False
    for l in open(pdbqt):
        if l.startswith("MODEL"):
            if inmodel:
                break
            inmodel = True
        elif l.startswith("ENDMDL"):
            break
        elif l.startswith("REMARK VINA RESULT") and aff is None:
            try:
                aff = float(l.split()[3])
            except (IndexError, ValueError):
                pass
        elif l.startswith(("ATOM", "HETATM")):
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            typ.append(l[77:79].strip() or l[12:16].strip()[:2])
    return np.array(xyz), np.array(typ, dtype=object), aff


def angle(a, b, c):
    v1, v2 = a - b, c - b
    cos = float(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def icc(groups):
    groups = [g for g in groups if len(g) >= 2]
    if len(groups) < 3:
        return None
    k = float(np.mean([len(g) for g in groups])); n = len(groups)
    grand = float(np.mean(np.concatenate([np.asarray(g) for g in groups])))
    msb = float(sum(len(g) * (np.mean(g) - grand) ** 2 for g in groups) / (n - 1))
    dfw = sum(len(g) - 1 for g in groups)
    msw = float(sum(((np.asarray(g) - np.mean(g)) ** 2).sum()
                    for g in groups) / dfw)
    s2b = max(0.0, (msb - msw) / k)
    return (s2b / (s2b + msw) if (s2b + msw) > 0 else 0.0,
            math.sqrt(s2b), math.sqrt(msw), n)


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    d = math.sqrt(float((rx ** 2).sum()) * float((ry ** 2).sum()))
    return float(rx @ ry / d) if d else float("nan")


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    x = x - x.mean(); y = y - y.mean()
    d = math.sqrt(float((x ** 2).sum()) * float((y ** 2).sum()))
    return float(x @ y / d) if d else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docked", required=True)
    ap.add_argument("--receptors", required=True)
    ap.add_argument("--prep", required=True)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--endpoints", default=None)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--end", type=int, default=537)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    # ---- alignment frame; also a FULL residue -> Torpedo position map ------
    aln = {}
    for h, s in read_fasta(a.aln):
        aln.setdefault(h.split()[0], s)
    tor = next((k for k in aln if "Torpedo" in k), None)
    res2col, r, col = {}, 0, 0
    for ch in aln[tor]:
        if ch == "-": col += 1
        elif ch == ".": continue
        elif ch.isupper(): r += 1; res2col[r] = col; col += 1
        else: r += 1; res2col[r] = None
    raw = "".join(c for c in aln[tor] if c not in "-.").upper()
    offset = (ELBOW.search(raw).start() + 3) - 200
    col2mature = {c: up - offset for up, c in res2col.items() if c is not None}

    def col_of(mat):
        up = mat + offset
        c = res2col.get(up)
        while c is None and up in res2col:
            up += 1; c = res2col.get(up)
        return c
    c_start, c_end = col_of(a.start), col_of(a.end)
    site_cols = {l: col_of(int(re.search(r"\d+", l).group())) for l in SITES}

    def mapping(sid):
        """trimmed residue index -> Torpedo mature number, plus site lookups."""
        idx2mat, sites, kept, cur = {}, {}, 0, 0
        want = {v: k for k, v in site_cols.items()}
        for ch in aln.get(sid, ""):
            if ch == ".": continue
            if ch == "-": cur += 1; continue
            if c_start <= cur <= c_end:
                kept += 1
                if ch.isupper():
                    if cur in col2mature:
                        idx2mat[kept] = col2mature[cur]
                    if cur in want:
                        sites[want[cur]] = kept
            if ch.isupper(): cur += 1
        return idx2mat, sites

    prep = {r["species"]: r for r in read_tsv(a.prep) if r.get("status") == "ok"}
    ends = {}
    for r in read_tsv(a.endpoints) if a.endpoints else []:
        sp = (r.get("Species_representative") or r.get("species") or
              r.get("Species_Name") or "").replace(" ", "_")
        y = r.get("log10_LC50", "")
        if sp and y not in ("", None):
            try:
                ends[sp] = float(y)
            except ValueError:
                pass
    print(f"{len(prep)} prepped receptors | {len(ends)} species with an endpoint")

    # ---- 1. affinities ----------------------------------------------------
    per_seed, best_pose = [], {}
    for spdir in sorted(glob.glob(os.path.join(a.docked, "*"))):
        if not os.path.isdir(spdir):
            continue
        sp = os.path.basename(spdir)
        for pq in sorted(glob.glob(os.path.join(spdir, "*_seed*.pdbqt"))):
            b = os.path.basename(pq)[:-6]
            lig, seed = b.rsplit("_seed", 1)
            xyz, typ, aff = first_pose(pq)
            if aff is None or not len(xyz):
                continue
            per_seed.append({"species": sp, "ligand": lig, "seed": int(seed),
                             "affinity": aff, "n_atoms": len(xyz)})
            key = (sp, lig)
            if key not in best_pose or aff < best_pose[key][0]:
                best_pose[key] = (aff, xyz, typ, pq)
    if not per_seed:
        sys.exit(f"no poses parsed under {a.docked}")

    with open(os.path.join(a.outdir, "affinity_per_seed.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["species", "ligand", "seed",
                                       "affinity", "n_atoms"])
        w.writeheader(); w.writerows(per_seed)

    by = defaultdict(list)
    for r in per_seed:
        by[(r["species"], r["ligand"])].append(r["affinity"])
    sp_aff = []
    for (sp, lig), v in sorted(by.items()):
        sp_aff.append({"species": sp, "ligand": lig, "n_seeds": len(v),
                       "affinity_mean": round(float(np.mean(v)), 3),
                       "affinity_sd": round(float(np.std(v, ddof=1)), 3)
                       if len(v) > 1 else 0.0,
                       "affinity_best": round(float(np.min(v)), 3)})
    with open(os.path.join(a.outdir, "affinity_species.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=list(sp_aff[0]))
        w.writeheader(); w.writerows(sp_aff)

    ligs = sorted({r["ligand"] for r in per_seed})
    var_rows = []
    print("\n== 1. is the affinity about species, or about the search seed? ==")
    print(f"{'ligand':26s} {'ICC':>6s} {'between_sd':>11s} {'within_sd':>10s} "
          f"{'n_sp':>5s}")
    for lig in ligs:
        groups = [v for (sp, l), v in by.items() if l == lig]
        res = icc(groups)
        if res is None:
            continue
        i, sb, sw, n = res
        var_rows.append({"ligand": lig, "ICC": round(i, 3),
                         "between_species_sd": round(sb, 3),
                         "within_species_sd": round(sw, 3), "n_species": n,
                         "usable": "yes" if i >= 0.75 else
                                   ("marginal" if i >= 0.5 else "NO_mostly_noise")})
        print(f"{lig:26s} {i:6.3f} {sb:11.3f} {sw:10.3f} {n:5d}")
    if var_rows:
        with open(os.path.join(a.outdir, "affinity_variance.tsv"), "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                               fieldnames=list(var_rows[0]))
            w.writeheader(); w.writerows(var_rows)
    else:
        print("  (fewer than 3 species with replicate seeds; ICC not estimable)")

    # ---- 2 + 3. geometry and contacts, on the best oxon pose --------------
    nac, contacts = [], defaultdict(Counter)
    for sp, p in sorted(prep.items()):
        key = (sp, OXON)
        if key not in best_pose:
            continue
        aff, lxyz, ltyp, pq = best_pose[key]
        rec = os.path.join(a.receptors, f"{sp}.pdb")
        if not os.path.exists(rec):
            continue
        pxyz, presi, presn, pname = pdb_heavy(rec)
        idx2mat, sites = mapping(p["sequence_id"])

        row = {"species": sp, "affinity": aff, "pose": os.path.basename(pq)}
        og = pxyz[(presi == int(p["ser_resnum"])) & (pname == "OG")]
        P = lxyz[ltyp == "P"]
        if len(og) and len(P):
            og, P = og[0], P[0]
            row["d_SerOG_P"] = round(float(np.linalg.norm(P - og)), 2)
            # oxygens bonded to P
            oidx = [i for i, t in enumerate(ltyp) if t in ("OA", "O")
                    and np.linalg.norm(lxyz[i] - P) < 1.85]
            arom = [i for i, t in enumerate(ltyp) if t == "A"]
            lg = None
            for i in oidx:
                if any(np.linalg.norm(lxyz[i] - lxyz[j]) < 1.55 for j in arom):
                    lg = i; break          # ester O joining P to the pyridine
            if lg is not None:
                row["angle_OG_P_Olg"] = round(angle(og, P, lxyz[lg]), 1)
                row["d_P_Olg"] = round(float(np.linalg.norm(lxyz[lg] - P)), 2)
                # deviation from the 180 deg in-line geometry
                row["inline_deviation"] = round(abs(180.0 -
                                                    row["angle_OG_P_Olg"]), 1)
            # oxyanion hole: backbone N of Gly118/Gly119/Ala201 to the P=O
            po = [i for i in oidx
                  if not any(np.linalg.norm(lxyz[i] - lxyz[j]) < 1.6
                             for j in range(len(lxyz)) if j != i
                             and not np.array_equal(lxyz[j], P))]
            oxy = []
            for lab in OXYANION:
                if lab in sites:
                    n_at = pxyz[(presi == sites[lab]) & (pname == "N")]
                    if len(n_at) and po:
                        oxy.append(min(float(np.linalg.norm(n_at[0] - lxyz[i]))
                                       for i in po))
            if oxy:
                row["oxyanion_min_dist"] = round(min(oxy), 2)
                row["n_oxyanion_hbond"] = sum(1 for d in oxy if d <= 3.5)
        # contacts, in Torpedo numbering
        if len(lxyz):
            d = np.linalg.norm(pxyz[:, None, :] - lxyz[None, :, :], axis=2)
            near = presi[(d <= CONTACT_CUT).any(axis=1)]
            mats = sorted({idx2mat[int(i)] for i in set(near.tolist())
                           if int(i) in idx2mat})
            row["n_contact_residues"] = len(mats)
            row["contacts_torpedo"] = ";".join(str(m) for m in mats)
            for m in mats:
                contacts[m][sp] = 1
        nac.append(row)

    if nac:
        cols = ["species", "affinity", "d_SerOG_P", "angle_OG_P_Olg",
                "inline_deviation", "d_P_Olg", "oxyanion_min_dist",
                "n_oxyanion_hbond", "n_contact_residues", "contacts_torpedo",
                "pose"]
        with open(os.path.join(a.outdir, "nac_geometry.tsv"), "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                               fieldnames=cols, extrasaction="ignore")
            w.writeheader(); w.writerows(nac)
        dd = [r["d_SerOG_P"] for r in nac if "d_SerOG_P" in r]
        an = [r["inline_deviation"] for r in nac if "inline_deviation" in r]
        print(f"\n== 2. near-attack geometry, {OXON}, {len(nac)} species ==")
        if dd:
            print(f"  Ser Og to P      : {min(dd):.2f} to {max(dd):.2f} A "
                  f"(median {sorted(dd)[len(dd)//2]:.2f})")
            n_close = sum(1 for x in dd if x <= 4.0)
            print(f"  within 4 A of P  : {n_close}/{len(dd)} species")
        if an:
            print(f"  deviation from in-line 180 deg: {min(an):.0f} to "
                  f"{max(an):.0f} (median {sorted(an)[len(an)//2]:.0f})")
            print("  A pose far from 4 A and far from in-line is not a "
                  "phosphorylation-competent complex, whatever its score.")

        n_sp = len(nac)
        crows = [{"torpedo_position": m, "n_species": len(v),
                  "frac_species": round(len(v) / n_sp, 3)}
                 for m, v in sorted(contacts.items(), key=lambda x: -len(x[1]))]
        if not crows:
            crows = [{"torpedo_position": "", "n_species": 0,
                      "frac_species": 0.0}]
        with open(os.path.join(a.outdir, "contact_frequency.tsv"), "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                               fieldnames=list(crows[0]))
            w.writeheader(); w.writerows(crows)
        print(f"\n== 3. contacts in Torpedo numbering ==")
        print("  most frequently contacted positions:")
        for c in crows[:12]:
            lab = next((s for s in SITES
                        if int(re.search(r"\d+", s).group()) == c["torpedo_position"]),
                       "")
            print(f"    {c['torpedo_position']:5d} {lab:8s} "
                  f"{c['n_species']:3d}/{n_sp} species")
        var = [c for c in crows if 0.15 <= c["frac_species"] <= 0.85]
        print(f"  positions contacted in some species but not others: {len(var)}")
        print("  Those are where the interaction differs across species, and "
              "they are the residue-level hypotheses worth testing against "
              "sensitivity.")

    # ---- 4. control diagnostic -------------------------------------------
    if ends:
        aff_by = defaultdict(dict)
        for r in sp_aff:
            aff_by[r["ligand"]][r["species"]] = r["affinity_mean"]
        print(f"\n== 4. CONTROL DIAGNOSTIC: affinity against log10 LC50 ==")
        crow = []
        for lig in ligs:
            pairs = [(aff_by[lig][s], ends[s]) for s in aff_by[lig] if s in ends]
            if len(pairs) < 6:
                continue
            x = [p[0] for p in pairs]; y = [p[1] for p in pairs]
            crow.append({"ligand": lig, "n_species": len(pairs),
                         "pearson_r": round(pearson(x, y), 3),
                         "spearman_rho": round(spearman(x, y), 3)})
        crow.sort(key=lambda r: -abs(r["spearman_rho"]))
        if not crow:
            print("  fewer than 6 species with both an affinity and an "
                  "endpoint; correlation not computed.")
        else:
            with open(os.path.join(a.outdir, "ligand_lc50_correlation.tsv"), "w",
                      newline="") as fh:
                w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                                   fieldnames=list(crow[0]))
                w.writeheader(); w.writerows(crow)
        if crow:
            print(f"{'ligand':26s} {'n':>4s} {'pearson':>8s} {'spearman':>9s}")
        for r in crow:
            print(f"{r['ligand']:26s} {r['n_species']:4d} "
                  f"{r['pearson_r']:8.3f} {r['spearman_rho']:9.3f}")
        d = {r["ligand"]: abs(r["spearman_rho"]) for r in crow}
        if OXON in d:
            for ctrl in ("chlorpyrifos", "TCP", "TCP_anion"):
                if ctrl in d and d[ctrl] >= d[OXON] - 0.05:
                    print(f"\n  WARNING: {ctrl} predicts LC50 as well as the "
                          f"oxon ({d[ctrl]:.2f} vs {d[OXON]:.2f}). {ctrl} "
                          "cannot phosphorylate acetylcholinesterase. The "
                          "ranking is therefore not about target engagement, "
                          "and that is a finding about the METHOD which must "
                          "be reported, not a nuisance to be dropped.")

    print(f"\n  {a.outdir}/")
    print("\nNothing here is a potency. Vina scores the non-covalent complex; "
          "OP potency is k2/Kd. Treat the geometry and the oxon-versus-control "
          "contrasts as the results, and carry the affinity ICC into any model "
          "that uses the scores.")


if __name__ == "__main__":
    main()
