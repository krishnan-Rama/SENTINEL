#!/usr/bin/env python3
"""
Receptor PDBQT conversion, with a check that the histidine work survived it.

The trap
--------
8B-prep_receptors.py established that all 82 catalytic triads are correctly oriented
and renamed the catalytic histidine HID so that Ne2 is free to accept the proton
from Ser Og. Conversion tools then add hydrogens according to their own rules.
Open Babel does not recognise HID as a histidine tautomer directive; asked to
protonate at pH 7.4 it will place the proton wherever its own logic says, and
can put it on Ne2, which blocks the acceptor and un-activates the nucleophile in
every structure. The careful work upstream is undone silently.

So this script converts, then INSPECTS the result: it looks for a polar hydrogen
on Nd1 and the absence of one on Ne2 for the catalytic histidine in each output,
and reports any structure where that is not the case. A receptor failing this
check should not be docked.

  module load micromamba/2.8.1
  export PATH=/mnt/ecotox/GROUP-smbpk/c23048124/envs/dock/bin:$PATH
  python3 8C-convert_receptors.py --prepped RECEPTORS/prepped \\
      --prep-table RECEPTORS/receptor_prep.tsv --outdir RECEPTORS/pdbqt
"""

import argparse, csv, glob, math, os, shutil, subprocess, sys


def read_atoms(path):
    out = []
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            try:
                out.append({"name": l[12:16].strip(), "resn": l[17:20].strip(),
                            "resi": int(l[22:26]),
                            "xyz": (float(l[30:38]), float(l[38:46]),
                                    float(l[46:54]))})
            except ValueError:
                continue
    return out


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepped", required=True)
    ap.add_argument("--prep-table", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ph", type=float, default=7.4)
    ap.add_argument("--rename-hid", action="store_true", default=True,
                    help="write HIS in the file handed to Open Babel, which "
                         "does not understand HID, then verify the tautomer "
                         "afterwards from the hydrogen positions")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    tmpd = os.path.join(a.outdir, "_tmp"); os.makedirs(tmpd, exist_ok=True)

    if not shutil.which("obabel"):
        sys.exit("obabel not on PATH. Run setup_dock_env.sh first.")

    prep = {}
    with open(a.prep_table, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("status") == "ok":
                prep[r["species"]] = r

    rows = []
    for pdb in sorted(glob.glob(os.path.join(a.prepped, "*.pdb"))):
        sp = os.path.basename(pdb)[:-4]
        if sp not in prep:
            print(f"  skip {sp}: not 'ok' in the prep table"); continue
        p = prep[sp]
        his_i = int(p["his_resnum"]); ser_i = int(p["ser_resnum"])
        row = {"species": sp, "his_resnum": his_i, "ser_resnum": ser_i}

        # Open Babel needs the standard residue name to treat it as histidine
        src = os.path.join(tmpd, f"{sp}.pdb")
        with open(pdb) as fi, open(src, "w") as fo:
            for l in fi:
                if l.startswith("ATOM") and l[17:20].strip() == "HID":
                    l = l[:17] + "HIS" + l[20:]
                fo.write(l)

        out = os.path.join(a.outdir, f"{sp}.pdbqt")
        r = subprocess.run(["obabel", src, "-opdbqt", "-O", out,
                            "-xr", "-p", str(a.ph)],
                           capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(out):
            row["status"] = "FAIL_conversion"
            row["detail"] = (r.stderr or "")[-200:]
            rows.append(row); print(f"  FAIL {sp}: conversion"); continue

        at = read_atoms(out)
        his = [x for x in at if x["resi"] == his_i]
        nd1 = next((x for x in his if x["name"] == "ND1"), None)
        ne2 = next((x for x in his if x["name"] == "NE2"), None)
        og = next((x for x in at if x["resi"] == ser_i and x["name"] == "OG"),
                  None)
        if not (nd1 and ne2 and og):
            row["status"] = "FAIL_triad_atoms_lost_in_conversion"
            rows.append(row); print(f"  FAIL {sp}: triad atoms lost"); continue

        # polar H are retained in PDBQT; assign each to its nearest ring N
        hs = [x for x in his if x["name"].startswith("H")]
        h_nd1 = sum(1 for h in hs if dist(h["xyz"], nd1["xyz"]) < 1.3)
        h_ne2 = sum(1 for h in hs if dist(h["xyz"], ne2["xyz"]) < 1.3)
        row.update(n_polarH_on_ND1=h_nd1, n_polarH_on_NE2=h_ne2,
                   d_NE2_SerOG=round(dist(ne2["xyz"], og["xyz"]), 2),
                   n_atoms=len(at))

        if h_ne2 > 0 and h_nd1 == 0:
            row["status"] = "FAIL_wrong_tautomer_HIE"
        elif h_ne2 > 0 and h_nd1 > 0:
            row["status"] = "FAIL_protonated_HIP"
        elif h_nd1 > 0:
            row["status"] = "ok"
        else:
            # no polar H written on either ring nitrogen
            row["status"] = "REVIEW_no_ring_H"
        rows.append(row)
        print(f"  {row['status']:26s} {sp:32s} "
              f"H(ND1)={h_nd1} H(NE2)={h_ne2} NE2-OG={row['d_NE2_SerOG']}")

    cols = ["species", "status", "his_resnum", "ser_resnum",
            "n_polarH_on_ND1", "n_polarH_on_NE2", "d_NE2_SerOG", "n_atoms",
            "detail"]
    with open(os.path.join(a.outdir, "receptor_pdbqt_qc.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols,
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    from collections import Counter
    print("\n== conversion status ==")
    for k, v in Counter(r["status"] for r in rows).most_common():
        print(f"  {k}: {v}")
    bad = [r for r in rows if r["status"] != "ok"]
    if bad:
        print(f"\n{len(bad)} receptors are NOT ready to dock.")
        print("FAIL_wrong_tautomer_HIE means Open Babel put the proton on Ne2, "
              "which blocks the acceptor that takes the proton from Ser Og. "
              "Docking into those measures a dead enzyme.")
        print("If most structures fail, switch the protonation step to pdb2pqr "
              "(it is in the env) which honours explicit tautomer assignment, "
              "and convert the resulting structure rather than the raw PDB.")
    print(f"\n  {a.outdir}/  ({sum(1 for r in rows if r['status']=='ok')} ready)"
          f"\n  {a.outdir}/receptor_pdbqt_qc.tsv")


if __name__ == "__main__":
    main()
