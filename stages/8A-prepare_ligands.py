#!/usr/bin/env python3
"""
Build the ligand panel, with identity verification.

Every SMILES here is checked against an independently known molecular formula
and monoisotopic-free molecular weight before anything is written. A wrong
SMILES that still parses would dock happily and produce a full set of plausible
numbers for the wrong molecule, and nothing downstream would catch it.

The panel and why each member is in it
--------------------------------------
  chlorpyrifos_oxon   The actual inhibitor. Everything else is a control.
  chlorpyrifos        NEGATIVE CONTROL. A phosphorothioate; it needs CYP450
                      desulfuration before it can phosphorylate anything. If it
                      predicts LC50 as well as the oxon does, the pipeline is
                      ranking lipophilicity, not target engagement, and that is
                      a diagnostic you want to have built in deliberately.
  TCP                 SECOND NEGATIVE CONTROL, the leaving group. Small,
                      chlorinated, aromatic, lipophilic, and chemically
                      incapable of phosphorylating a serine. If species rank
                      similarly on TCP and on the oxon, your ranking is about
                      the gorge accommodating a chlorinated aromatic rather
                      than about organophosphate chemistry. This is the control
                      I would least like to see omitted.
  TCP_anion           TCP has a pKa near 4.5, so above 99% is anionic at pH 7.4.
                      Both forms are built because the neutral is the membrane
                      permeant species and the anion is what is present in
                      solution; docking both bounds the answer.
  acetylcholine       The natural substrate, permanently cationic. This is the
                      competitive reference: the BARMIE logic of "can the
                      chemical outcompete the endogenous ligand" only makes
                      sense with the endogenous ligand in the comparison.
  desethyl_CPO        A CPO metabolite, anionic phosphate monoester.

  module load micromamba/2.8.1
  export PATH=/mnt/ecotox/GROUP-smbpk/c23048124/envs/dock/bin:$PATH
  python3 8A-prepare_ligands.py --outdir LIGANDS
"""

import argparse, csv, os, re, subprocess, sys

# name: (SMILES, expected formula, expected MW, net charge, role)
PANEL = {
    "chlorpyrifos_oxon": (
        "CCOP(=O)(OCC)Oc1nc(Cl)c(Cl)cc1Cl",
        "C9H11Cl3NO4P", 334.52, 0, "inhibitor"),
    "chlorpyrifos": (
        "CCOP(=S)(OCC)Oc1nc(Cl)c(Cl)cc1Cl",
        "C9H11Cl3NO3PS", 350.59, 0, "negative_control_parent"),
    "TCP": (
        "Oc1nc(Cl)c(Cl)cc1Cl",
        "C5H2Cl3NO", 198.43, 0, "negative_control_leaving_group"),
    "TCP_anion": (
        "[O-]c1nc(Cl)c(Cl)cc1Cl",
        "C5HCl3NO", 197.43, -1, "negative_control_leaving_group_pH7.4"),
    "acetylcholine": (
        "CC(=O)OCC[N+](C)(C)C",
        "C7H16NO2", 146.21, +1, "natural_substrate_reference"),
    "desethyl_CPO": (
        "CCOP(=O)(O)Oc1nc(Cl)c(Cl)cc1Cl",
        "C7H7Cl3NO4P", 306.47, 0, "metabolite"),
}
MW_TOL = 0.6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--conformers", type=int, default=20)
    ap.add_argument("--mw-tol", type=float, default=MW_TOL)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        sys.exit("rdkit not found. Run setup_dock_env.sh and put its bin on PATH.")

    rows, ok_names = [], []
    for name, (smi, formula, mw, chg, role) in PANEL.items():
        row = {"ligand": name, "smiles": smi, "role": role,
               "expected_formula": formula, "expected_mw": mw,
               "expected_charge": chg}
        m = Chem.MolFromSmiles(smi)
        if m is None:
            row["status"] = "FAIL_smiles_unparseable"
            rows.append(row); print(f"  FAIL {name}: SMILES will not parse"); continue

        got_f = rdMolDescriptors.CalcMolFormula(m)
        got_mw = Descriptors.MolWt(m)
        got_c = Chem.GetFormalCharge(m)
        row.update(got_formula=got_f, got_mw=round(got_mw, 2), got_charge=got_c)

        # RDKit appends a charge suffix, e.g. C7H16NO2+ or C5HCl3NO-. Strip
        # only that suffix. rstrip("+-0123456789") is wrong: on C7H16NO2+ it
        # removes the '+' and then keeps going and removes the '2', silently
        # turning a correct formula into a mismatch.
        fbase = re.sub(r"[+-]\d*$", "", got_f)
        problems = []
        if fbase != formula:
            problems.append(f"formula {got_f} != {formula}")
        if abs(got_mw - mw) > a.mw_tol:
            problems.append(f"MW {got_mw:.2f} != {mw}")
        if got_c != chg:
            problems.append(f"charge {got_c} != {chg}")
        if problems:
            row["status"] = "FAIL_identity: " + "; ".join(problems)
            rows.append(row)
            print(f"  FAIL {name}: {'; '.join(problems)}")
            continue

        mh = Chem.AddHs(m)
        ps = AllChem.ETKDGv3(); ps.randomSeed = 0xC0FFEE
        cids = AllChem.EmbedMultipleConfs(mh, numConfs=a.conformers, params=ps)
        if not len(cids):
            row["status"] = "FAIL_embed"; rows.append(row)
            print(f"  FAIL {name}: 3D embedding failed"); continue
        res = AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=2000)
        energies = [e for conv, e in res]
        best = int(min(range(len(energies)), key=lambda i: energies[i]))

        sdf = os.path.join(a.outdir, f"{name}.sdf")
        w = Chem.SDWriter(sdf); w.write(mh, confId=cids[best]); w.close()

        pdbqt = os.path.join(a.outdir, f"{name}.pdbqt")
        r = subprocess.run(["mk_prepare_ligand.py", "-i", sdf, "-o", pdbqt],
                           capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(pdbqt):
            r2 = subprocess.run(["obabel", sdf, "-opdbqt", "-O", pdbqt],
                                capture_output=True, text=True)
            row["pdbqt_via"] = "obabel" if r2.returncode == 0 else "FAILED"
        else:
            row["pdbqt_via"] = "meeko"

        nrot = sum(1 for l in open(pdbqt) if l.startswith("BRANCH")) \
            if os.path.exists(pdbqt) else 0
        row.update(n_conformers=len(cids), mmff_energy=round(energies[best], 2),
                   n_rotatable_branches=nrot, sdf=sdf, pdbqt=pdbqt,
                   status="ok" if row.get("pdbqt_via") != "FAILED"
                          else "FAIL_pdbqt")
        rows.append(row)
        ok_names.append(name)
        print(f"  ok   {name:20s} {got_f:16s} MW {got_mw:7.2f} "
              f"q{got_c:+d}  {nrot} torsions  via {row['pdbqt_via']}")

    cols = ["ligand", "role", "status", "smiles", "expected_formula",
            "got_formula", "expected_mw", "got_mw", "expected_charge",
            "got_charge", "n_conformers", "mmff_energy",
            "n_rotatable_branches", "pdbqt_via", "sdf", "pdbqt"]
    with open(os.path.join(a.outdir, "ligand_manifest.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols,
                           extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    print(f"\n{len(ok_names)}/{len(PANEL)} ligands built -> {a.outdir}")
    bad = [r for r in rows if r["status"] != "ok"]
    if bad:
        print("NOT BUILT:")
        for r in bad:
            print(f"  {r['ligand']}: {r['status']}")
        print("\nA failed identity check means the SMILES is wrong, not that "
              "the molecule is difficult. Fix it before docking; a wrong "
              "structure produces a complete and entirely plausible set of "
              "results.")
    print(f"  {a.outdir}/ligand_manifest.tsv")


if __name__ == "__main__":
    main()
