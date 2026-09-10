#!/usr/bin/env python3
"""
Receptor preparation: catalytic triad orientation and the docking box.

The two things that silently ruin comparative OP docking
--------------------------------------------------------
1. THE CATALYTIC HISTIDINE IS OFTEN FLIPPED. In a serine hydrolase triad,
   His Ne2 accepts the proton from Ser Og, and His Nd1 donates to the Glu
   carboxylate. Structure prediction (and crystallography) frequently place the
   imidazole ring 180 degrees out, because Nd1/Ne2 and Cd2/Ce1 are nearly
   isosteric and the density or the model cannot tell them apart. If Ne2 is not
   the nitrogen facing Ser Og, the nucleophile is not activated in the structure
   you dock into, and every score in the study inherits the error.

   This script measures both distances, flips the ring where they are reversed,
   and renames the residue HID (neutral, proton on Nd1) so downstream tools do
   not reassign it.

2. THE BOX IS USUALLY CENTRED ON THE WRONG THING. Centring on the protein, or
   on the catalytic serine alone, either misses the gorge or clips the
   trichloropyridinol leaving group, which must point out toward the gorge
   mouth for a phosphorylation-competent pose. Here the box is built along the
   actual gorge axis: from the Ser Og at the base to the centroid of the
   peripheral anionic site (Trp279, Tyr70, Asp72, Tyr121) at the mouth,
   computed per structure. Gorge geometry differs between species; a single
   shared box would impose one species' geometry on all of them.

Outputs a corrected PDB per species, a per-structure box, and a report you must
read before docking.

  module purge
  module load Python/3.11.5-GCCcore-13.2.0 SciPy-bundle/2025.07-gfbf-2025b
  python3 8B-prep_receptors.py --models QC/best_models --aln ../CLASSIFY/aligned.a2m \\
      --qc QC/model_qc.tsv --outdir RECEPTORS
"""

import argparse, csv, glob, math, os, re, sys
from collections import defaultdict

try:
    import numpy as np
except ImportError:
    sys.exit("numpy required. module load SciPy-bundle/2025.07-gfbf-2025b")

ELBOW = re.compile(r"G[EQ]S[AG]G")
SITES = ["Tyr70", "Asp72", "Trp84", "Gly118", "Gly119", "Tyr121", "Glu199",
         "Ser200", "Ala201", "Trp279", "Phe288", "Phe290", "Glu327",
         "Phe330", "Tyr334", "His440"]
PAS = ["Trp279", "Tyr70", "Asp72", "Tyr121"]      # gorge mouth
# Imidazole atom pairs exchanged by a 180 degree ring flip about CB-CG.
HIS_FLIP = [("ND1", "CD2"), ("CE1", "NE2")]
HBOND_MAX = 3.6                                    # generous; models are not xtal


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
    """list of dicts, preserving line order and formatting fields."""
    at = []
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            at.append({"line": l, "serial": l[6:11], "name": l[12:16].strip(),
                       "alt": l[16], "resn": l[17:20].strip(), "chain": l[21],
                       "resi": int(l[22:26]),
                       "xyz": np.array([float(l[30:38]), float(l[38:46]),
                                        float(l[46:54])])})
    return at


def write_pdb(path, atoms):
    with open(path, "w") as f:
        for a in atoms:
            l = a["line"]
            x, y, z = a["xyz"]
            f.write(f"{l[:17]}{a['resn']:>3}{l[20:30]}"
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}")
        f.write("END\n")


def coord(atoms, resi, name):
    for a in atoms:
        if a["resi"] == resi and a["name"] == name:
            return a["xyz"]
    return None


def sidechain_centroid(atoms, resi):
    pts = [a["xyz"] for a in atoms if a["resi"] == resi
           and a["name"] not in ("N", "CA", "C", "O")]
    return np.mean(pts, axis=0) if pts else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", required=True)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--qc", default=None)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--end", type=int, default=537)
    ap.add_argument("--pad", type=float, default=6.0,
                    help="padding beyond the gorge axis endpoints, Angstrom")
    ap.add_argument("--min-box", type=float, default=20.0)
    a = ap.parse_args()
    prep = os.path.join(a.outdir, "prepped")
    os.makedirs(prep, exist_ok=True)

    # ---- alignment frame, identical to trimming and QC --------------------
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
    site_cols = {l: col_of(int(re.search(r"\d+", l).group())) for l in SITES}
    print(f"offset {offset}; trimmed columns {c_start}-{c_end}")

    # sequence_id per species, from the QC table if available
    seqid = {}
    if a.qc and os.path.exists(a.qc):
        with open(a.qc, newline="") as fh:
            for rr in csv.DictReader(fh, delimiter="\t"):
                seqid[rr["species"]] = rr["sequence_id"]

    def site_positions(sid):
        if sid not in aln:
            return {}
        out, kept, cur = {}, 0, 0
        want = {v: k for k, v in site_cols.items()}
        for ch in aln[sid]:
            if ch == ".":
                continue
            if ch == "-":
                cur += 1; continue
            if c_start <= cur <= c_end:
                kept += 1
                if ch.isupper() and cur in want:
                    out[want[cur]] = kept
            if ch.isupper():
                cur += 1
        return out

    rows = []
    for pdb in sorted(glob.glob(os.path.join(a.models, "*.pdb"))):
        sp = os.path.basename(pdb)[:-4]
        sid = seqid.get(sp)
        pos = site_positions(sid) if sid else {}
        row = {"species": sp, "sequence_id": sid or ""}
        if not all(k in pos for k in ("Ser200", "His440", "Glu327")):
            row["status"] = "FAIL_triad_not_mapped"
            rows.append(row); continue

        at = read_pdb(pdb)
        s_i, h_i, e_i = pos["Ser200"], pos["His440"], pos["Glu327"]
        row.update(ser_resnum=s_i, his_resnum=h_i, glu_resnum=e_i)

        og = coord(at, s_i, "OG")
        nd1, ne2 = coord(at, h_i, "ND1"), coord(at, h_i, "NE2")
        glu_o = [c for c in (coord(at, e_i, "OE1"), coord(at, e_i, "OE2"))
                 if c is not None]
        resn_h = next((x["resn"] for x in at if x["resi"] == h_i), "?")
        resn_s = next((x["resn"] for x in at if x["resi"] == s_i), "?")
        row["his_resname_in"] = resn_h
        if og is None or nd1 is None or ne2 is None or not glu_o:
            row["status"] = "FAIL_missing_triad_atoms"
            rows.append(row); continue
        if resn_s != "SER":
            row["status"] = f"FAIL_res200_is_{resn_s}"
            rows.append(row); continue

        d = lambda p, q: float(np.linalg.norm(p - q))
        dmin = lambda p, L: min(d(p, q) for q in L)
        # correct orientation: NE2 to Ser Og, ND1 to Glu carboxylate
        as_is = d(ne2, og) + dmin(nd1, glu_o)
        flipped = d(nd1, og) + dmin(ne2, glu_o)
        do_flip = flipped < as_is

        if do_flip:
            for n1, n2 in HIS_FLIP:
                A = [x for x in at if x["resi"] == h_i and x["name"] == n1]
                B = [x for x in at if x["resi"] == h_i and x["name"] == n2]
                if A and B:
                    A[0]["xyz"], B[0]["xyz"] = B[0]["xyz"].copy(), A[0]["xyz"].copy()
            nd1, ne2 = coord(at, h_i, "ND1"), coord(at, h_i, "NE2")

        row["his_flipped"] = bool(do_flip)
        row["d_NE2_SerOG"] = round(d(ne2, og), 2)
        row["d_ND1_GluO"] = round(dmin(nd1, glu_o), 2)
        row["triad_ok"] = bool(row["d_NE2_SerOG"] <= HBOND_MAX
                               and row["d_ND1_GluO"] <= HBOND_MAX)

        # neutral, proton on Nd1, so Ne2 can accept from Ser Og
        for x in at:
            if x["resi"] == h_i and x["resn"] in ("HIS", "HIE", "HIP", "HSD",
                                                  "HSE", "HSP"):
                x["resn"] = "HID"
        row["his_resname_out"] = "HID"

        # ---- gorge axis and box ------------------------------------------
        mouth_pts = [sidechain_centroid(at, pos[p]) for p in PAS if p in pos]
        mouth_pts = [p for p in mouth_pts if p is not None]
        row["n_pas_residues"] = len(mouth_pts)
        if not mouth_pts:
            row["status"] = "FAIL_no_PAS_residues"
            rows.append(row); continue
        mouth = np.mean(mouth_pts, axis=0)
        axis = mouth - og
        gorge_len = float(np.linalg.norm(axis))
        centre = og + axis * 0.5
        # box must contain both endpoints plus room for the leaving group
        span = gorge_len + 2 * a.pad
        size = max(span, a.min_box)
        row.update(gorge_length=round(gorge_len, 2),
                   box_center_x=round(float(centre[0]), 3),
                   box_center_y=round(float(centre[1]), 3),
                   box_center_z=round(float(centre[2]), 3),
                   box_size_x=round(size, 1), box_size_y=round(size, 1),
                   box_size_z=round(size, 1),
                   ser_og_x=round(float(og[0]), 3),
                   ser_og_y=round(float(og[1]), 3),
                   ser_og_z=round(float(og[2]), 3))
        row["status"] = "ok" if row["triad_ok"] else "REVIEW_triad_geometry"

        write_pdb(os.path.join(prep, f"{sp}.pdb"), at)
        rows.append(row)

    cols = ["species", "sequence_id", "status", "ser_resnum", "his_resnum",
            "glu_resnum", "his_resname_in", "his_resname_out", "his_flipped",
            "d_NE2_SerOG", "d_ND1_GluO", "triad_ok", "n_pas_residues",
            "gorge_length", "ser_og_x", "ser_og_y", "ser_og_z",
            "box_center_x", "box_center_y", "box_center_z",
            "box_size_x", "box_size_y", "box_size_z"]
    with open(os.path.join(a.outdir, "receptor_prep.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols,
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    from collections import Counter
    print("\n== status ==")
    for k, v in Counter(r.get("status", "?") for r in rows).most_common():
        print(f"  {k}: {v}")
    fl = [r for r in rows if r.get("his_flipped")]
    print(f"\ncatalytic histidines flipped: {len(fl)}/{len(rows)}")
    if fl:
        print("  " + ", ".join(r["species"] for r in fl[:12])
              + (" ..." if len(fl) > 12 else ""))
        print("  Each of these would have docked with an unactivated "
              "nucleophile had the ring not been corrected.")
    ok = [r for r in rows if r.get("status") == "ok"]
    if ok:
        gl = [r["gorge_length"] for r in ok]
        dn = [r["d_NE2_SerOG"] for r in ok]
        print(f"\nprepared: {len(ok)}")
        print(f"  Ne2 to Ser Og : {min(dn):.2f} to {max(dn):.2f} A "
              f"(median {sorted(dn)[len(dn)//2]:.2f})")
        print(f"  gorge length  : {min(gl):.1f} to {max(gl):.1f} A "
              f"(median {sorted(gl)[len(gl)//2]:.1f})")
        print("\n  Gorge length is your first cross-species structural "
              "measurement. Compare its spread against the all-atom seed RMSD "
              "from QC before treating any of it as a species difference.")
    bad = [r for r in rows if r.get("status") != "ok"]
    for r in bad:
        print(f"  {r['species']:32s} {r.get('status')} "
              f"NE2-OG={r.get('d_NE2_SerOG','')} ND1-GluO={r.get('d_ND1_GluO','')}")
    print(f"\n  {prep}/\n  {a.outdir}/receptor_prep.tsv")


if __name__ == "__main__":
    main()
