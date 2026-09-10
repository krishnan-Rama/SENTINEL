#!/usr/bin/env python3
"""
Bi-phasic analysis of the AChE gorge, from poses you already have.

The mechanism
-------------
The AChE active site is a narrow gorge with two ligand binding sites: a
peripheral site (P-site) at the mouth and an acylation site (A-site) at the
base. Organophosphates are captured at the P-site, descend the gorge, and
phosphorylate the catalytic serine at the A-site. Between them sits a
free-energy basin at Trp84 (W86 in mammalian numbering) formed by cation-pi
interactions, which is where an unconstrained docking search will settle
because it is the thermodynamic minimum of the NON-COVALENT complex.

That is exactly what your first analysis found: median Ser Og to P of 7.9 A,
Trp84 contacted in 81 of 82 species, and no pose within 4 A of the serine. The
docking did not fail; it reported the basin the literature says is there. The
reactive pose is a higher-energy transient state and a rigid docking score will
never rank it first.

What this script does
---------------------
Rather than taking the single best pose, it reads EVERY mode from EVERY seed
(up to 27 per species-ligand pair), projects each onto the gorge axis, and
scores the phases separately:

  dG_P      best affinity among poses at the gorge mouth
  dG_M      best affinity in the mid-gorge basin
  dG_A      best affinity at the acylation site
  ddG_PA    dG_P minus dG_A, a per-species proxy for whether the gorge favours
            capture or descent

ddG_PA is mechanistically closer to k_i than any single affinity, because for
chlorpyrifos-oxon it is Kd, not k2, that varies, and the two-site occupancy is
what makes Kd vary.

No new docking is required.

  module purge
  module load Python/3.11.5-GCCcore-13.2.0 SciPy-bundle/2025.07-gfbf-2025b
  python3 9B-analyse_biphasic.py --docked DOCKED --prep RECEPTORS/receptor_prep.tsv \\
      --endpoints ../REPORT/analysis_set.tsv --outdir BIPHASIC
"""

import argparse, csv, glob, math, os, re, sys
from collections import defaultdict, Counter

try:
    import numpy as np
except ImportError:
    sys.exit("numpy required. module load SciPy-bundle/2025.07-gfbf-2025b")
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_PLT = True
except ImportError:
    HAVE_PLT = False

OXON = "chlorpyrifos_oxon"
A_FRAC, P_FRAC = 0.33, 0.66      # legacy fractional cuts, retained for comparison

# THREE CORRECTIONS to the first version of this analysis, each of which
# changed the headline result:
#
# 1. ANCHOR ATOM. The original anchored a pose at the phosphorus when the
#    ligand had one and at the whole-molecule centroid otherwise. Phosphorus
#    sits at one END of chlorpyrifos-oxon; a centroid sits in the MIDDLE. The
#    apparent difference in gorge penetration between acetylcholine and the
#    oxon split perfectly along which anchor was used, so the comparison was
#    between two reference points rather than between two molecules. All
#    ligands are now anchored identically. Default is the minimum distance from
#    ANY ligand heavy atom to the catalytic serine, which answers "did any part
#    of the ligand reach the acylation site" for every chemistry equally.
#
# 2. ZONE BOUNDARIES. Cuts at fixed fractions of gorge length give a species
#    with a longer gorge a physically larger acylation zone, so the penetration
#    fraction partly restates gorge length (observed r up to +0.76). Cuts are
#    now in absolute Angstrom from the catalytic serine by default.
#
# 3. POSE COUNTS ARE NOT OCCUPANCY. Vina returns diversity-filtered local
#    minima under an RMSD cutoff, not Boltzmann-weighted samples, and mode rank
#    predicts zone. A count of poses per zone measures the search, not the
#    thermodynamics. Counts are still written but renamed, and the energy-based
#    quantities are what should be interpreted.
A_ABS, P_ABS = 6.0, 10.0         # A-site <= 6 A from Ser Og; P-site >= 10 A


def read_tsv(p, required=True):
    if not os.path.exists(p):
        if required:
            sys.exit(f"missing: {p}")
        return []
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def all_poses(pdbqt):
    """Every MODEL in a Vina output: (affinity, coords, autodock types)."""
    out, aff, xyz, typ = [], None, [], []
    for l in open(pdbqt):
        if l.startswith("MODEL"):
            aff, xyz, typ = None, [], []
        elif l.startswith("REMARK VINA RESULT") and aff is None:
            try:
                aff = float(l.split()[3])
            except (IndexError, ValueError):
                pass
        elif l.startswith(("ATOM", "HETATM")):
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            typ.append(l[77:79].strip() or l[12:16].strip()[:2])
        elif l.startswith("ENDMDL"):
            if aff is not None and xyz:
                out.append((aff, np.array(xyz), np.array(typ, dtype=object)))
            aff, xyz, typ = None, [], []
    if aff is not None and xyz:
        out.append((aff, np.array(xyz), np.array(typ, dtype=object)))
    return out


def icc(groups):
    """Between-species vs within-species (between-seed) variance, as ICC."""
    groups = [g for g in groups if len(g) >= 2]
    if len(groups) < 3:
        return None
    k = float(np.mean([len(g) for g in groups])); n = len(groups)
    grand = float(np.mean(np.concatenate([np.asarray(g, float) for g in groups])))
    msb = float(sum(len(g) * (np.mean(g) - grand) ** 2 for g in groups) / (n - 1))
    dfw = sum(len(g) - 1 for g in groups)
    if dfw <= 0:
        return None
    msw = float(sum(((np.asarray(g, float) - np.mean(g)) ** 2).sum()
                    for g in groups) / dfw)
    s2b = max(0.0, (msb - msw) / k)
    return ((s2b / (s2b + msw) if (s2b + msw) > 0 else 0.0),
            math.sqrt(s2b), math.sqrt(msw), n)


def angle(a, b, c):
    v1, v2 = a - b, c - b
    n = np.linalg.norm(v1) * np.linalg.norm(v2)
    if n == 0:
        return float("nan")
    return math.degrees(math.acos(max(-1.0, min(1.0, float(v1 @ v2 / n)))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--docked", required=True)
    ap.add_argument("--prep", required=True)
    ap.add_argument("--endpoints", default=None)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--anchor", default="min",
                    choices=["min", "centroid", "phosphorus"],
                    help="reference point for pose depth, applied identically "
                         "to every ligand. 'min' = closest ligand heavy atom "
                         "to the catalytic serine.")
    ap.add_argument("--zones", default="absolute",
                    choices=["absolute", "fractional"],
                    help="absolute uses --a-cut/--p-cut in Angstrom; "
                         "fractional reproduces the original, flawed, "
                         "gorge-length-dependent behaviour")
    ap.add_argument("--a-cut", type=float, default=A_ABS)
    ap.add_argument("--p-cut", type=float, default=P_ABS)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    print(f"anchor = {a.anchor} (identical for all ligands); "
          f"zones = {a.zones}"
          + (f" (A <= {a.a_cut} A, P >= {a.p_cut} A)" if a.zones == "absolute"
             else " (fractions of gorge length)"))

    prep = {}
    for r in read_tsv(a.prep):
        if r.get("status") != "ok":
            continue
        og = np.array([fnum(r["ser_og_x"]), fnum(r["ser_og_y"]),
                       fnum(r["ser_og_z"])])
        ctr = np.array([fnum(r["box_center_x"]), fnum(r["box_center_y"]),
                        fnum(r["box_center_z"])])
        # the box centre is the midpoint of Ser Og -> PAS centroid, so the
        # mouth is recovered by doubling the vector
        mouth = og + 2.0 * (ctr - og)
        L = float(np.linalg.norm(mouth - og))
        if L < 1e-6:
            continue
        prep[r["species"]] = {"og": og, "mouth": mouth, "L": L,
                              "u": (mouth - og) / L,
                              "gorge_length": fnum(r.get("gorge_length")) or L}
    print(f"{len(prep)} receptors with a gorge axis")

    ends = {}
    for r in read_tsv(a.endpoints, required=False) if a.endpoints else []:
        sid = r.get("Sequence_id", "")
        tag = sid.split("__")[0] if "__" in sid else (
            r.get("Species_representative", "").replace(" ", "_"))
        v = fnum(r.get("log10_LC50"))
        if tag and v is not None:
            ends[tag] = {"y": v, "class": r.get("Class", ""),
                         "name": r.get("Species_representative", tag)}
    print(f"{len(ends)} species with an endpoint")

    # ---- 1. every pose, projected onto the gorge axis ---------------------
    poses = []
    for spdir in sorted(glob.glob(os.path.join(a.docked, "*"))):
        sp = os.path.basename(spdir)
        if not os.path.isdir(spdir) or sp not in prep:
            continue
        g = prep[sp]
        for pq in sorted(glob.glob(os.path.join(spdir, "*_seed*.pdbqt"))):
            b = os.path.basename(pq)[:-6]
            lig, seed = b.rsplit("_seed", 1)
            for mode, (aff, xyz, typ) in enumerate(all_poses(pq), 1):
                P = xyz[typ == "P"]
                # identical anchoring for every ligand; see notes at the top
                if a.anchor == "centroid":
                    anchor = xyz.mean(axis=0)
                elif a.anchor == "phosphorus" and len(P):
                    anchor = P[0]
                else:
                    anchor = xyz[int(np.argmin(
                        np.linalg.norm(xyz - g["og"], axis=1)))]
                t = float((anchor - g["og"]) @ g["u"])
                frac = t / g["L"]
                d_og = float(np.linalg.norm(anchor - g["og"]))
                if a.zones == "absolute":
                    zone = ("A" if d_og <= a.a_cut else
                            ("P" if d_og >= a.p_cut else "M"))
                else:
                    zone = ("A" if frac < A_FRAC else
                            ("P" if frac > P_FRAC else "M"))
                rec = {"species": sp, "ligand": lig, "seed": int(seed),
                       "mode": mode, "affinity": aff,
                       "depth": round(t, 2), "frac_depth": round(frac, 4),
                       "zone": zone, "has_P": bool(len(P)),
                       "anchor_mode": a.anchor,
                       "d_SerOG_anchor": round(d_og, 2),
                       # only meaningful where a phosphorus exists; the old
                       # column name implied it always did
                       "d_SerOG_P": (round(float(np.linalg.norm(P[0] - g["og"])), 2)
                                     if len(P) else "")}
                if len(P):
                    Pp = P[0]
                    oidx = [i for i, ty in enumerate(typ) if ty in ("OA", "O")
                            and np.linalg.norm(xyz[i] - Pp) < 1.85]
                    arom = [i for i, ty in enumerate(typ) if ty == "A"]
                    lg = next((i for i in oidx
                               if any(np.linalg.norm(xyz[i] - xyz[j]) < 1.55
                                      for j in arom)), None)
                    if lg is not None:
                        rec["inline_deviation"] = round(
                            abs(180.0 - angle(g["og"], Pp, xyz[lg])), 1)
                poses.append(rec)
    if not poses:
        sys.exit(f"no poses parsed under {a.docked}")
    print(f"{len(poses)} poses parsed "
          f"({len({(p['species'], p['ligand']) for p in poses})} pairs)")

    cols = ["species", "ligand", "seed", "mode", "affinity", "depth",
            "frac_depth", "zone", "d_SerOG_P", "inline_deviation", "has_P"]
    with open(os.path.join(a.outdir, "poses_all.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(poses)

    print("\n== pose occupancy by zone ==")
    for lig in sorted({p["ligand"] for p in poses}):
        sub = [p for p in poses if p["ligand"] == lig]
        c = Counter(p["zone"] for p in sub)
        n = len(sub)
        print(f"  {lig:22s} A {c['A']*100//n:3d}%  M {c['M']*100//n:3d}%  "
              f"P {c['P']*100//n:3d}%   (n={n})")
    print("  A = acylation site, M = mid-gorge basin, P = peripheral site.")
    print("  NOTE: these are counts of Vina output modes, which are "
          "diversity-filtered local minima, not a Boltzmann ensemble. Treat "
          "them as a description of the search, and interpret the per-zone "
          "ENERGIES below rather than the counts.")

    # ---- 2. best affinity per zone ---------------------------------------
    best = defaultdict(dict)
    for p in poses:
        k = (p["species"], p["ligand"])
        z = p["zone"]
        if z not in best[k] or p["affinity"] < best[k][z]["affinity"]:
            best[k][z] = p
    rows = []
    for (sp, lig), zs in sorted(best.items()):
        r = {"species": sp, "ligand": lig}
        for z in ("A", "M", "P"):
            if z in zs:
                r[f"dG_{z}"] = zs[z]["affinity"]
                r[f"n_{z}"] = sum(1 for p in poses if p["species"] == sp
                                  and p["ligand"] == lig and p["zone"] == z)
        if "dG_P" in r and "dG_A" in r:
            r["ddG_PA"] = round(r["dG_P"] - r["dG_A"], 3)
        if "A" in zs:
            r["A_d_SerOG_P"] = zs["A"].get("d_SerOG_P", "")
            r["A_inline_deviation"] = zs["A"].get("inline_deviation", "")
        if sp in ends:
            r["log10_LC50"] = ends[sp]["y"]; r["class"] = ends[sp]["class"]
        rows.append(r)
    fields = ["species", "ligand", "class", "log10_LC50", "dG_A", "dG_M",
              "dG_P", "ddG_PA", "n_A", "n_M", "n_P", "A_d_SerOG_P",
              "A_inline_deviation"]
    with open(os.path.join(a.outdir, "biphasic_species.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    ox = [r for r in rows if r["ligand"] == OXON]
    have = [r for r in ox if "ddG_PA" in r]
    print(f"\n== {OXON}: phases resolved for {len(have)}/{len(ox)} species ==")
    if have:
        for k, lab in (("dG_A", "acylation site"), ("dG_M", "mid-gorge basin"),
                       ("dG_P", "peripheral site")):
            v = [r[k] for r in have if k in r]
            if v:
                print(f"  {lab:18s} {min(v):6.2f} to {max(v):6.2f} "
                      f"(median {sorted(v)[len(v)//2]:.2f}) kcal/mol")
        d = [r["ddG_PA"] for r in have]
        print(f"  ddG_PA (P minus A) {min(d):+.2f} to {max(d):+.2f} "
              f"(median {sorted(d)[len(d)//2]:+.2f})")
        print("  Positive means the acylation site is the better-scoring of "
              "the two; negative means the mouth is.")

    # ---- 2b. penetration descriptors, with their own noise floor ---------
    # Each seed contributes 9 modes, so the per-seed A-site fraction is a
    # replicate measurement of the same species property. That gives an ICC on
    # the same footing as the structural descriptors: a penetration fraction
    # whose between-seed variance rivals its between-species variance is not a
    # species trait and must not be regressed on.
    per_seed = defaultdict(list)
    for pp in poses:
        per_seed[(pp["species"], pp["ligand"], pp["seed"])].append(pp["zone"])
    seedfrac = defaultdict(lambda: defaultdict(list))
    for (sp, lig, sd), zs in per_seed.items():
        if zs:
            seedfrac[lig][sp].append(sum(1 for z in zs if z == "A") / len(zs))

    ligs_all = sorted(seedfrac)
    pen_rows, icc_rows = defaultdict(dict), []
    for lig in ligs_all:
        res = icc(list(seedfrac[lig].values()))
        if res:
            i, sb, sw, nsp = res
            icc_rows.append({"descriptor": f"pen_modefracA_{lig}", "n_species": nsp,
                             "between_species_sd": round(sb, 4),
                             "within_species_sd": round(sw, 4),
                             "ICC": round(i, 3),
                             "usable": "yes" if i >= 0.75 else
                                       ("marginal" if i >= 0.5 else "NO_mostly_noise")})
        for sp, v in seedfrac[lig].items():
            # renamed: it is the fraction of returned MODES, not an occupancy
            pen_rows[sp][f"pen_modefracA_{lig}"] = round(float(np.mean(v)), 4)
            pen_rows[sp][f"pen_modefracA_{lig}_sd"] = round(
                float(np.std(v, ddof=1)) if len(v) > 1 else 0.0, 4)
    for r in rows:
        if "ddG_PA" in r:
            pen_rows[r["species"]][f"pen_ddGPA_{r['ligand']}"] = r["ddG_PA"]

    keys = sorted({k for d in pen_rows.values() for k in d})
    with open(os.path.join(a.outdir, "penetration_descriptors.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["species"] + keys, extrasaction="ignore")
        w.writeheader()
        for sp in sorted(pen_rows):
            w.writerow({"species": sp, **pen_rows[sp]})
    if icc_rows:
        with open(os.path.join(a.outdir, "penetration_icc.tsv"), "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                               fieldnames=list(icc_rows[0]))
            w.writeheader(); w.writerows(icc_rows)
        print("\n== penetration fraction: species trait or seed noise? ==")
        for r in sorted(icc_rows, key=lambda x: -x["ICC"]):
            print(f"  {r['descriptor']:34s} ICC {r['ICC']:5.3f}  "
                  f"between {r['between_species_sd']:.3f}  "
                  f"within {r['within_species_sd']:.3f}   {r['usable']}")

    # ---- 3. figures -------------------------------------------------------
    if not HAVE_PLT:
        print("\nmatplotlib unavailable; figures skipped")
        print(f"\n  {a.outdir}/"); return

    classes = sorted({e["class"] for e in ends.values() if e["class"]})
    pal = plt.get_cmap("tab10")
    ccol = {c: pal(i % 10) for i, c in enumerate(classes)}

    # fig 1: docking energy landscape along the gorge axis
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axes[0]
    for lig, col in zip([OXON, "TCP", "acetylcholine", "chlorpyrifos"],
                        ["#c0392b", "#2c7fb8", "#27ae60", "#8e44ad"]):
        sub = [p for p in poses if p["ligand"] == lig]
        if not sub:
            continue
        f = np.array([p["frac_depth"] for p in sub])
        e = np.array([p["affinity"] for p in sub])
        bins = np.linspace(-0.1, 1.2, 22)
        cen, mn = [], []
        for i in range(len(bins) - 1):
            m = (f >= bins[i]) & (f < bins[i + 1])
            if m.sum() >= 5:
                cen.append(0.5 * (bins[i] + bins[i + 1]))
                mn.append(float(np.percentile(e[m], 10)))
        if cen:
            ax.plot(cen, mn, "-o", ms=3, lw=1.4, color=col, label=lig)
    ax.axvspan(-0.1, A_FRAC, color="0.85", alpha=.5)
    ax.axvspan(P_FRAC, 1.2, color="0.85", alpha=.5)
    ax.text(0.10, ax.get_ylim()[1], "A-site", fontsize=8, va="top")
    ax.text(0.90, ax.get_ylim()[1], "P-site", fontsize=8, va="top")
    ax.set_xlabel("fractional depth along gorge axis (0 = Ser O$\\gamma$)")
    ax.set_ylabel("best decile affinity (kcal/mol)")
    ax.set_title("A  docking energy profile along the gorge", fontsize=10,
                 loc="left")
    ax.legend(fontsize=7, frameon=False)

    ax = axes[1]
    for lig in sorted({p["ligand"] for p in poses}):
        sub = [p["frac_depth"] for p in poses if p["ligand"] == lig]
        ax.hist(sub, bins=np.linspace(-0.1, 1.2, 30), histtype="step", lw=1.3,
                label=lig, density=True)
    ax.axvline(A_FRAC, ls=":", c="0.5"); ax.axvline(P_FRAC, ls=":", c="0.5")
    ax.set_xlabel("fractional depth"); ax.set_ylabel("pose density")
    ax.set_title("B  where the poses actually sit", fontsize=10, loc="left")
    ax.legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "fig1_gorge_landscape.png"), dpi=200)
    fig.savefig(os.path.join(a.outdir, "fig1_gorge_landscape.svg"))
    plt.close(fig)

    # fig 2: P-site versus A-site affinity, per species
    if have:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
        ax = axes[0]
        for r in have:
            ax.scatter(r["dG_A"], r["dG_P"], s=34, edgecolors="0.3", lw=.4,
                       color=ccol.get(r.get("class", ""), "0.6"))
        lim = [min(min(r["dG_A"] for r in have), min(r["dG_P"] for r in have)) - .3,
               max(max(r["dG_A"] for r in have), max(r["dG_P"] for r in have)) + .3]
        ax.plot(lim, lim, ls="--", c="0.6", lw=.9)
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("$\\Delta$G acylation site (kcal/mol)")
        ax.set_ylabel("$\\Delta$G peripheral site (kcal/mol)")
        ax.set_title("A  two-site affinity, chlorpyrifos-oxon", fontsize=10,
                     loc="left")
        ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=ccol[c],
                                      label=c, markersize=6) for c in classes],
                  fontsize=7, frameon=False, loc="upper left")

        ax = axes[1]
        wy = [r for r in have if "log10_LC50" in r]
        if wy:
            for r in wy:
                ax.scatter(r["ddG_PA"], r["log10_LC50"], s=38,
                           edgecolors="0.3", lw=.4,
                           color=ccol.get(r.get("class", ""), "0.6"))
            x = np.array([r["ddG_PA"] for r in wy])
            yv = np.array([r["log10_LC50"] for r in wy])
            if len(x) > 3 and x.std() > 0:
                b, c0 = np.polyfit(x, yv, 1)
                xs = np.linspace(x.min(), x.max(), 50)
                ax.plot(xs, b * xs + c0, c="0.4", lw=1.1)
                rr = float(np.corrcoef(x, yv)[0, 1])
                ax.set_title(f"B  $\\Delta\\Delta$G$_{{P-A}}$ vs sensitivity "
                             f"(r = {rr:+.2f}, n = {len(x)})",
                             fontsize=10, loc="left")
        ax.set_xlabel("$\\Delta\\Delta$G$_{P-A}$ (kcal/mol)")
        ax.set_ylabel("log10 LC50 (mg/L)")
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "fig2_two_site.png"), dpi=200)
        fig.savefig(os.path.join(a.outdir, "fig2_two_site.svg"))
        plt.close(fig)

    # fig 3: near-attack geometry of the best A-site poses
    asite = [p for p in poses if p["zone"] == "A" and p["ligand"] == OXON
             and p.get("d_SerOG_P") not in ("", None)
             and "inline_deviation" in p]
    if asite:
        fig, ax = plt.subplots(figsize=(6.2, 5))
        d = np.array([float(p["d_SerOG_P"]) for p in asite])
        ang = np.array([float(p["inline_deviation"]) for p in asite])
        e = np.array([p["affinity"] for p in asite])
        s = ax.scatter(d, ang, c=e, cmap="viridis_r", s=26, edgecolors="none")
        ax.axvline(4.0, ls="--", c="#c0392b", lw=1.1)
        ax.axhline(30.0, ls="--", c="#c0392b", lw=1.1)
        ax.add_patch(plt.Rectangle((0, 0), 4.0, 30.0, fc="#c0392b", alpha=.10))
        ax.text(0.3, 34, "near-attack conformation", fontsize=8,
                color="#c0392b")
        ax.set_xlabel("Ser O$\\gamma$ to P (${\\rm \\AA}$)")
        ax.set_ylabel("deviation from in-line 180$^\\circ$ ($^\\circ$)")
        ax.set_title(f"A-site poses, {OXON} ({len(asite)} poses)", fontsize=10)
        fig.colorbar(s, ax=ax, label="affinity (kcal/mol)")
        n_nac = int(((d <= 4.0) & (ang <= 30.0)).sum())
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "fig3_near_attack.png"), dpi=200)
        fig.savefig(os.path.join(a.outdir, "fig3_near_attack.svg"))
        plt.close(fig)
        print(f"\n  A-site poses in a near-attack conformation "
              f"(<=4 A and <=30 deg): {n_nac}/{len(asite)}")
        if n_nac == 0:
            print("  Even restricted to the acylation site, no pose reaches "
                  "reactive geometry. That is the quantitative case for "
                  "covalent docking or MD, and it is a result rather than a "
                  "gap: a rigid non-covalent search cannot sample the "
                  "transition-state-like arrangement.")

    print(f"\n  {a.outdir}/poses_all.tsv"
          f"\n  {a.outdir}/biphasic_species.tsv"
          f"\n  fig1_gorge_landscape, fig2_two_site, fig3_near_attack "
          "(.png and .svg)")
    print(f"\n  {a.outdir}/penetration_descriptors.tsv"
          f"\n  {a.outdir}/penetration_icc.tsv")
    print("\nFeed these to 10B-model_sensitivity.py:")
    print(f"  --penetration {a.outdir}/penetration_descriptors.tsv \\")
    print(f"  --penetration-icc {a.outdir}/penetration_icc.tsv")
    print("\nAdd dG_A, dG_P and ddG_PA to 10B-model_sensitivity.py as predictors. "
          "ddG_PA is the one worth watching: for chlorpyrifos-oxon it is Kd "
          "rather than k2 that varies with concentration, and two-site "
          "occupancy is what makes Kd vary.")


if __name__ == "__main__":
    main()
