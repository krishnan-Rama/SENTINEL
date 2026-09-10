#!/usr/bin/env python3
"""
Structural descriptors of the AChE gorge, with an honest noise model.

Why this computes everything five times
---------------------------------------
Your gorge lengths span 12.2 to 14.5 A across 82 species, SD 0.66 A. Your
between-seed all-atom RMSD at the active site runs from 0.22 to 1.50 A. For the
noisiest species the prediction noise is comparable to the entire between-species
range. A descriptor computed from one model per species cannot distinguish a
species difference from a seed difference, and no amount of downstream
statistics repairs that.

So every descriptor is computed on all five relaxed models per species, and the
variance is partitioned: sigma^2 between species against sigma^2 within species
(between seeds). The intraclass correlation, ICC = s2b / (s2b + s2w), is the
fraction of the descriptor's variance that is actually about species.

  ICC > 0.75   usable; species differences dominate prediction noise
  ICC 0.5-0.75 marginal; report with the ICC and expect wide CIs in the model
  ICC < 0.5    the descriptor is mostly seed noise. Do NOT put it in a PGLS
               and then interpret its coefficient.

Reporting the ICC alongside every structural predictor is the thing that
separates this from a comparative docking paper that computed one model each
and never checked.

Descriptors
-----------
  gorge_length        Ser Og to the peripheral anionic site centroid
  gorge_volume        probe-accessible volume in the gorge channel, grid based
  bottleneck_radius   narrowest inscribed radius along the gorge axis
  mouth_radius        inscribed radius at the gorge mouth
  n_aromatic_lining   aromatic residues lining the channel
  net_charge_lining   formal charge of the lining, relevant to the
                      electrostatic steering of the cationic substrate
  acyl_pocket_width   Phe288 to Phe290 side-chain centroid separation
  trp84_trp279_dist   choline site to peripheral site span

  module purge
  module load Python/3.11.5-GCCcore-13.2.0 SciPy-bundle/2025.07-gfbf-2025b
  python3 7-structural_descriptors.py --models FOLD/models \\
      --aln ../CLASSIFY/aligned.a2m --qc QC/model_qc.tsv --outdir STRUCT
"""

import argparse, csv, glob, math, os, re, sys
from collections import defaultdict, Counter

try:
    import numpy as np
    from scipy.spatial import cKDTree
    from scipy import ndimage
except ImportError:
    sys.exit("numpy+scipy required. module load SciPy-bundle/2025.07-gfbf-2025b")

ELBOW = re.compile(r"G[EQ]S[AG]G")
SITES = ["Tyr70", "Asp72", "Trp84", "Gly118", "Gly119", "Tyr121", "Glu199",
         "Ser200", "Ala201", "Trp279", "Phe288", "Phe290", "Glu327",
         "Phe330", "Tyr334", "His440"]
PAS = ["Trp279", "Tyr70", "Asp72", "Tyr121"]
AROMATIC = {"PHE", "TYR", "TRP", "HIS", "HID", "HIE", "HIP"}
POS = {"ARG", "LYS"}
NEG = {"ASP", "GLU"}
PROBE = 3.0          # heavy-atom centre to grid point; ~vdW 1.6 + probe 1.4
GRID = 0.7           # A
CYL_R = 9.0          # channel radius considered
LINING_R = 8.0       # residue heavy atom within this of the axis


def read_fasta(p):
    recs, h, s = [], None, []
    for l in open(p):
        if l.startswith(">"):
            if h: recs.append((h, "".join(s)))
            h, s = l[1:].rstrip("\n"), []
        elif l.strip(): s.append(l.strip())
    if h: recs.append((h, "".join(s)))
    return recs


def read_pdb(path):
    xyz, resi, resn, name = [], [], [], []
    for l in open(path):
        if l.startswith("ATOM"):
            nm = l[12:16].strip()
            if nm.startswith("H") or l[76:78].strip() == "H":
                continue
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            resi.append(int(l[22:26])); resn.append(l[17:20].strip()); name.append(nm)
    return (np.array(xyz), np.array(resi), np.array(resn, dtype=object),
            np.array(name, dtype=object))


def sc_centroid(xyz, resi, name, i):
    m = (resi == i) & ~np.isin(name, ["N", "CA", "C", "O"])
    return xyz[m].mean(0) if m.any() else None


def channel_descriptors(xyz, base, mouth):
    """Grid the cylinder between base and mouth; return volume and radii."""
    axis = mouth - base
    L = float(np.linalg.norm(axis))
    if L < 1e-6:
        return {}
    u = axis / L
    tree = cKDTree(xyz)

    # local frame
    tmp = np.array([1.0, 0, 0])
    if abs(u @ tmp) > 0.9:
        tmp = np.array([0, 1.0, 0])
    e1 = np.cross(u, tmp); e1 /= np.linalg.norm(e1)
    e2 = np.cross(u, e1)

    nt = int(L / GRID) + 1
    nr = int(CYL_R / GRID) + 1
    ts = np.linspace(0, L, nt)
    grid = np.zeros((nt, 2 * nr + 1, 2 * nr + 1), dtype=bool)
    coords = np.empty((nt, 2 * nr + 1, 2 * nr + 1, 3))
    offs = (np.arange(2 * nr + 1) - nr) * GRID
    for i, t in enumerate(ts):
        c = base + u * t
        A, B = np.meshgrid(offs, offs, indexing="ij")
        pts = c + A[..., None] * e1 + B[..., None] * e2
        coords[i] = pts
        inside = (A ** 2 + B ** 2) <= CYL_R ** 2
        d, _ = tree.query(pts.reshape(-1, 3))
        grid[i] = inside & (d.reshape(A.shape) > PROBE)

    # keep only the cavity connected to the axis near the base
    lab, n = ndimage.label(grid)
    if n == 0:
        return {"gorge_volume": 0.0}
    # Seed on the axis in the MIDDLE of the channel. Seeding near the base
    # fails systematically: the axis origin is the Ser Og atom itself, so every
    # voxel within one probe radius of the base is occupied by the very atom
    # that defines the axis, and the flood fill starts on a blocked voxel.
    # Collect every component the axis passes through between 25% and 75% and
    # take the largest, so a locally pinched channel does not zero the volume.
    lo, hi = int(0.25 * nt), max(int(0.75 * nt), int(0.25 * nt) + 1)
    axis_labs = [int(lab[i, nr, nr]) for i in range(lo, hi)
                 if lab[i, nr, nr] > 0]
    if not axis_labs:
        # axis fully blocked: fall back to the largest component anywhere in
        # the mid-section, and flag it by reporting the volume as usual
        mids = lab[lo:hi]
        cand = mids[mids > 0]
        if cand.size == 0:
            return {"gorge_volume": 0.0, "bottleneck_radius": 0.0,
                    "mouth_radius": 0.0, "base_radius": 0.0}
        axis_labs = [Counter(cand.tolist()).most_common(1)[0][0]]
    sizes = {l: int((lab == l).sum()) for l in set(axis_labs)}
    seed_lab = max(sizes, key=sizes.get)
    cav = lab == seed_lab
    vol = float(cav.sum()) * GRID ** 3

    # inscribed radius per slab: largest distance-to-protein among cavity pts
    radii = []
    for i in range(nt):
        m = cav[i]
        if not m.any():
            radii.append(0.0); continue
        d, _ = tree.query(coords[i][m])
        radii.append(float(d.max() - 1.4))     # subtract probe to get free radius
    radii = np.array(radii)
    mid = radii[int(0.15 * nt):int(0.85 * nt)]
    return {"gorge_volume": round(vol, 1),
            "bottleneck_radius": round(float(mid.min()), 3) if mid.size else 0.0,
            "mouth_radius": round(float(radii[-3:].mean()), 3),
            "base_radius": round(float(radii[:3].mean()), 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", required=True)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--qc", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--end", type=int, default=537)
    ap.add_argument("--max-models", type=int, default=5)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

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

    def col_of(mat):
        up = mat + offset
        c = res2col.get(up)
        while c is None and up in res2col:
            up += 1; c = res2col.get(up)
        return c
    c_start, c_end = col_of(a.start), col_of(a.end)
    site_cols = {l: col_of(int(re.search(r"\d+", l).group())) for l in SITES}

    def site_positions(sid):
        out, kept, cur = {}, 0, 0
        want = {v: k for k, v in site_cols.items()}
        for ch in aln.get(sid, ""):
            if ch == ".": continue
            if ch == "-": cur += 1; continue
            if c_start <= cur <= c_end:
                kept += 1
                if ch.isupper() and cur in want:
                    out[want[cur]] = kept
            if ch.isupper(): cur += 1
        return out

    keep = {}
    with open(a.qc, newline="") as fh:
        for rr in csv.DictReader(fh, delimiter="\t"):
            if rr["verdict"] == "pass":
                keep[rr["species"]] = rr["sequence_id"]
    print(f"{len(keep)} species passing QC")

    per_model = []
    for sp, sid in sorted(keep.items()):
        pos = site_positions(sid)
        if not all(k in pos for k in ("Ser200",)) or not any(p in pos for p in PAS):
            print(f"  skip {sp}: site mapping incomplete"); continue
        pdbs = sorted(glob.glob(os.path.join(
            a.models, "**", f"{sid}_relaxed_rank_00*_*.pdb"), recursive=True))
        for k, pdb in enumerate(pdbs[:a.max_models], 1):
            xyz, resi, resn, name = read_pdb(pdb)
            og = xyz[(resi == pos["Ser200"]) & (name == "OG")]
            if not len(og):
                continue
            og = og[0]
            mp = [sc_centroid(xyz, resi, name, pos[p]) for p in PAS if p in pos]
            mp = [x for x in mp if x is not None]
            if not mp:
                continue
            mouth = np.mean(mp, axis=0)
            row = {"species": sp, "sequence_id": sid, "model": k,
                   "gorge_length": round(float(np.linalg.norm(mouth - og)), 3)}
            row.update(channel_descriptors(xyz, og, mouth))

            # lining composition
            axis = mouth - og
            L = np.linalg.norm(axis); u = axis / L
            v = xyz - og
            t = np.clip(v @ u, 0, L)
            perp = np.linalg.norm(v - t[:, None] * u, axis=1)
            near = perp <= LINING_R
            lin = sorted(set(resi[near].tolist()))
            rn = {int(i): resn[resi == i][0] for i in lin}
            row["n_lining"] = len(lin)
            row["n_aromatic_lining"] = sum(1 for i in lin if rn[i] in AROMATIC)
            row["net_charge_lining"] = (sum(1 for i in lin if rn[i] in POS)
                                        - sum(1 for i in lin if rn[i] in NEG))
            row["frac_aromatic_lining"] = round(
                row["n_aromatic_lining"] / max(1, len(lin)), 4)

            if "Phe288" in pos and "Phe290" in pos:
                c1 = sc_centroid(xyz, resi, name, pos["Phe288"])
                c2 = sc_centroid(xyz, resi, name, pos["Phe290"])
                if c1 is not None and c2 is not None:
                    row["acyl_pocket_width"] = round(
                        float(np.linalg.norm(c1 - c2)), 3)
            if "Trp84" in pos and "Trp279" in pos:
                c1 = sc_centroid(xyz, resi, name, pos["Trp84"])
                c2 = sc_centroid(xyz, resi, name, pos["Trp279"])
                if c1 is not None and c2 is not None:
                    row["trp84_trp279_dist"] = round(
                        float(np.linalg.norm(c1 - c2)), 3)
            per_model.append(row)
        print(f"  {sp}: {len([r for r in per_model if r['species']==sp])} models")

    if not per_model:
        sys.exit("no descriptors computed")
    DESC = [c for c in ("gorge_length", "gorge_volume", "bottleneck_radius",
                        "mouth_radius", "base_radius", "n_lining",
                        "n_aromatic_lining", "frac_aromatic_lining",
                        "net_charge_lining", "acyl_pocket_width",
                        "trp84_trp279_dist")
            if any(c in r for r in per_model)]
    cols = ["species", "sequence_id", "model"] + DESC
    with open(os.path.join(a.outdir, "descriptors_per_model.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols,
                           extrasaction="ignore")
        w.writeheader(); w.writerows(per_model)

    # ---- species means and variance partition -----------------------------
    by_sp = defaultdict(list)
    for r in per_model:
        by_sp[r["species"]].append(r)

    sp_rows = []
    for sp, rs in sorted(by_sp.items()):
        row = {"species": sp, "n_models": len(rs)}
        for d in DESC:
            v = [r[d] for r in rs if d in r and r[d] is not None]
            if v:
                row[d] = round(float(np.mean(v)), 4)
                row[d + "_sd"] = round(float(np.std(v, ddof=1)), 4) if len(v) > 1 else 0.0
        sp_rows.append(row)
    scols = ["species", "n_models"] + [c for d in DESC for c in (d, d + "_sd")]
    with open(os.path.join(a.outdir, "descriptors_species.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=scols,
                           extrasaction="ignore")
        w.writeheader(); w.writerows(sp_rows)

    var_rows = []
    for d in DESC:
        groups = [[r[d] for r in rs if d in r] for rs in by_sp.values()]
        groups = [g for g in groups if len(g) >= 2]
        if len(groups) < 3:
            continue
        k = float(np.mean([len(g) for g in groups]))
        means = np.array([np.mean(g) for g in groups])
        grand = float(np.mean(np.concatenate(groups)))
        n = len(groups)
        msb = float(sum(len(g) * (np.mean(g) - grand) ** 2
                        for g in groups) / (n - 1))
        dfw = sum(len(g) - 1 for g in groups)
        msw = float(sum(((np.array(g) - np.mean(g)) ** 2).sum()
                        for g in groups) / dfw)
        s2b = max(0.0, (msb - msw) / k)
        icc = s2b / (s2b + msw) if (s2b + msw) > 0 else 0.0
        var_rows.append({
            "descriptor": d, "n_species": n,
            "between_species_sd": round(math.sqrt(s2b), 4),
            "within_species_sd": round(math.sqrt(msw), 4),
            "ICC": round(icc, 3),
            "species_range": f"{means.min():.3g} to {means.max():.3g}",
            "usable": "yes" if icc >= 0.75 else
                      ("marginal" if icc >= 0.5 else "NO_mostly_noise")})
    var_rows.sort(key=lambda r: -r["ICC"])
    with open(os.path.join(a.outdir, "variance_partition.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(var_rows[0]))
        w.writeheader(); w.writerows(var_rows)

    print("\n== variance partition: is the descriptor about species, or seeds? ==")
    print(f"{'descriptor':22s} {'ICC':>6s} {'between_sd':>11s} "
          f"{'within_sd':>10s}  species range")
    for r in var_rows:
        print(f"{r['descriptor']:22s} {r['ICC']:6.3f} "
              f"{r['between_species_sd']:11.3f} {r['within_species_sd']:10.3f}"
              f"  {r['species_range']}   {r['usable']}")

    good = [r["descriptor"] for r in var_rows if r["ICC"] >= 0.75]
    marg = [r["descriptor"] for r in var_rows if 0.5 <= r["ICC"] < 0.75]
    bad = [r["descriptor"] for r in var_rows if r["ICC"] < 0.5]
    print(f"\nusable (ICC >= 0.75): {', '.join(good) or 'NONE'}")
    print(f"marginal (0.5-0.75):  {', '.join(marg) or 'none'}")
    print(f"mostly seed noise:    {', '.join(bad) or 'none'}")
    if not good:
        print("\nNo descriptor clears 0.75. That is a real and reportable "
              "result: at AlphaFold-ensemble resolution the gorge geometry of "
              "these orthologues is not distinguishable between species by "
              "these measures. Do not proceed to interpret small structural "
              "differences, and say so.")
    print(f"\n  {a.outdir}/descriptors_per_model.tsv"
          f"\n  {a.outdir}/descriptors_species.tsv"
          f"\n  {a.outdir}/variance_partition.tsv")
    print("\nCarry the ICC into the modelling: a descriptor with ICC 0.6 has "
          "roughly 40% of its variance as noise, which attenuates its "
          "regression coefficient toward zero. Report the ICC beside every "
          "structural predictor rather than only its p-value.")


if __name__ == "__main__":
    main()
