# Stages

Every stage reads and writes tab-separated tables under `--outdir`. Nothing is held in
memory between stages, so any intermediate can be inspected, edited or replaced and the
downstream stages will pick up the change without anything upstream being re-run.

Paths below are relative to `--outdir`. The command shown is what `bin/sentinel` runs; you can
run it yourself with the same arguments.

---

## 1 anchors

Fetches the reference backbone by UniProt accession and verifies each entry against a length
window and a required substring of the protein name. A mismatch is a hard failure.

```
python3 stages/1-anchors.py --outdir S1_sequences --contact you@example.ac.uk
```

**Writes** `S1_sequences/anchors/reference_backbone.faa`, `reference_backbone.tsv`

The TSV carries a verdict column per anchor. If any row is not `pass`, stop. An unverified
backbone makes every downstream clade assignment unfalsifiable, and the failure mode is
silent: a wrong sequence produces a complete and plausible tree.

---

## 2 retrieve

Pulls the cholinesterase-family candidates from OrthoDB and builds the species inventory and
the ECOTOX coverage map.

**Reads** the cleaned ECOTOX table
**Writes** `ODB/by_species/`, `ODB/species_inventory.tsv`, `ODB/ecotox_coverage.tsv`

Retrieval is deliberately recall-oriented and does not filter on protein-name annotation.
In the reference dataset 41 per cent of entries had no curated description, and some
reference targets are named only by locus tag, so a text filter loses genuine orthologues
while keeping mislabelled ones.

*Script not in this release. See [below](#scripts-not-yet-released).*

---

## 3A classify

Aligns every candidate to a profile built from the verified anchors, establishes the residue
numbering frame by locating the nucleophile elbow motif, and scores 16 diagnostic
active-site positions.

```
python3 stages/3A-classify_active_site.py \
    --candidates ODB/by_species \
    --anchors S1_sequences/anchors/reference_backbone.faa \
    --outdir CLASSIFY --threads 16
```

**Writes** `CLASSIFY/active_site_calls.tsv`, `CLASSIFY/aligned.a2m`, `CLASSIFY/candidates.faa`

The numbering offset is derived, not assumed, and then validated: the script checks that
human AChE returns Phe/Phe at the acyl-pocket positions and human BChE returns Leu/Val, and
stops if it does not. Do not remove that check.

The `call` column is a convenience. The per-residue columns are the data, and the call is
superseded by the tree in stage 3B.

---

## 3B place

Builds one tree from every candidate plus the anchors, assigns each sequence to its nearest
anchor by patristic distance, reports the margin to the runner-up, and selects one
representative per species.

```
python3 stages/3B-assign_by_tree.py --classify CLASSIFY \
    --anchors S1_sequences/anchors/reference_backbone.faa \
    --outdir TREE --threads 16
```

**Writes** `TREE/tree_assignments.tsv`, `TREE/representatives.tsv`,
`TREE/representatives.faa`, `TREE/ache.tree`

Where the tree and the residue call disagree, both are reported and the row is flagged.
Resolve those by looking at them, not by dropping a column.

Known weakness, stated in the script: the anchor set contains only one non-insect
invertebrate reference, so molluscan, crustacean and annelid placements are weaker than
vertebrate or insect ones and should not be presented at equal confidence.

For a publication-quality tree with bootstraps, run the SLURM script instead of relying on
this one:

```
ROOT=$OUTDIR sbatch stages/3C-rebuild_tree.slurm
```

which realigns with L-INS-i, trims columns gapped in more than 60 per cent of sequences, and
runs RAxML-NG under LG+G4+F with 200 bootstraps into `PHYLO/ache.raxml.support`.

---

## 4 dataset

Joins the representatives back to the endpoints and splits the species into the two sets
that everything downstream keeps separate.

```
python3 stages/4A-finalise_endpoints.py       --ecotox $ECOTOX_TABLE --odb ODB --outdir FINAL
python3 stages/4B-build_final_set.py --final FINAL --tree TREE --outdir FINAL_SET
python3 stages/4C-build_master_table.py --classify CLASSIFY --tree TREE --odb ODB \
    --final FINAL --final-set FINAL_SET --ecotox $ECOTOX_TABLE --out master_table.csv
```

**Writes** `FINAL_SET/endpoint_final.tsv` and `.faa`, `FINAL_SET/prediction_final.tsv` and
`.faa`, `master_table.csv`

The endpoint set has a measured LC50 and a verified target: this is your n, and it is the
only set on which anything can be validated. The prediction set has a target and no
endpoint. Report them separately and never pool them into one performance figure.

Confidence tiers A, B and C are attached as a covariate rather than used as a filter. Tier C
means a judgement call was needed. Fit with and without.

*`4A-finalise_endpoints.py` is not in this release.*

---

## 5 fold

Trims every sequence to the catalytic domain at the same structural position, then predicts
five models per sequence.

```
python3 stages/5A-prepare_for_folding.py --aln CLASSIFY/aligned.a2m \
    --fasta FINAL_SET/endpoint_final.tsv FINAL_SET/prediction_final.tsv \
    --outdir FOLD --chunk 12
bin/sentinel stage 5B   # or: sbatch --array=0-N stages/5B-fold.slurm FOLD/chunks FOLD/models
```

**Writes** `FOLD/trim_report.tsv`, `FOLD/chunks/`, `FOLD/models/`

Trimming is not cosmetic. Folding the full precursor puts an unconstrained signal-peptide
helix somewhere, and in a fraction of models it lies across the gorge mouth, so any volume
or box computed from those models measures an artefact. It also drags whole-model pLDDT
around in a way that differs between species, which destroys the comparison you need.

Boundaries come from the shared profile alignment, not from a length heuristic, so every
sequence is cut at the same structural position.

Pass the TSVs, not the FASTAs. The FASTA headers carry only a species tag, and matching on
species tag selects every family member for that species rather than the chosen
representative. On the reference dataset that is roughly 650 sequences instead of 85. The
script refuses by default and tells you this.

Five models per sequence, six recycles, templates enabled, Amber relaxation. Templates
matter here: cholinesterase is exceptionally well covered in the PDB. `colabfold_batch`
skips sequences whose outputs already exist, so a timed-out array task resumes rather than
restarting.

---

## 6 qc

Measures the two things that decide whether a model is usable for this study.

```
python3 stages/6-qc_models.py --models FOLD/models --aln CLASSIFY/aligned.a2m \
    --final-set FINAL_SET --outdir QC
```

**Writes** `QC/model_qc.tsv`, `QC/best_models/`

Per-residue pLDDT at the 16 active-site positions, not the whole-model mean. A model can
average 92 with a ragged acyl pocket, and it will still dock and still return a number.

Between-seed active-site RMSD, backbone and all-atom. **Use the all-atom figure.** Gorge
volume and near-attack geometry are set by side-chain rotamers, not backbone, so the CA
number understates the uncertainty on exactly the quantities being compared across species.
Any between-species difference smaller than this is seed variance.

`ser_og_*` is the catalytic serine Og coordinate per model, and it is the docking box
centre. It is per-model. Do not reuse one species' coordinates for another.

---

## 7 descriptors

Computes gorge descriptors on all five models per species and partitions their variance.

```
python3 stages/7-structural_descriptors.py --models FOLD/models \
    --aln CLASSIFY/aligned.a2m --qc QC/model_qc.tsv --outdir STRUCT
```

**Writes** `STRUCT/descriptors_per_model.tsv`, `STRUCT/descriptors_species.tsv`,
`STRUCT/variance_partition.tsv`

`variance_partition.tsv` carries the intraclass correlation for each descriptor:

| ICC | Meaning |
|-----|---------|
| > 0.75 | usable; species differences dominate prediction noise |
| 0.5 to 0.75 | marginal; report the ICC and expect wide intervals |
| < 0.5 | mostly seed noise. Do not put it in a regression and interpret the coefficient |

---

## 8A ligands

Builds the compound panel from SMILES and verifies every molecule against an independently
known formula, weight and formal charge before writing anything.

```
python3 stages/8A-prepare_ligands.py --outdir LIGANDS
```

**Writes** `LIGANDS/*.sdf`, `LIGANDS/*.pdbqt`, `LIGANDS/ligand_manifest.tsv`

A wrong SMILES that still parses docks happily and produces a complete, plausible set of
numbers for the wrong molecule, and nothing downstream would catch it.

The panel must contain inactive analogues. In the reference panel these are the unactivated
parent compound, both protonation states of the leaving group, and the natural substrate.
They are the specificity control, and they are the reason the pipeline can falsify its own
premise.

---

## 8B, 8C receptors

Two steps, because two different things go wrong.

```
python3 stages/8B-prep_receptors.py --models QC/best_models --aln CLASSIFY/aligned.a2m \
    --qc QC/model_qc.tsv --outdir RECEPTORS
python3 stages/8C-convert_receptors.py --prepped RECEPTORS/prepped \
    --prep-table RECEPTORS/receptor_prep.tsv --outdir RECEPTORS/pdbqt
```

**Writes** `RECEPTORS/prepped/`, `RECEPTORS/receptor_prep.tsv`, `RECEPTORS/pdbqt/`,
`RECEPTORS/pdbqt/receptor_pdbqt_qc.tsv`

*Catalytic histidine orientation.* Nd1 and Ne2 are nearly isosteric, so prediction and
crystallography both place the imidazole ring 180 degrees out often enough that it must be
checked. If Ne2 does not face the catalytic serine, the nucleophile is not activated in the
structure you dock into, and every score inherits the error. The first script measures both
distances, flips the ring where they are reversed, and renames the residue HID.

*The box.* Centring on the protein misses the gorge; centring on the serine clips the
leaving group, which has to point toward the mouth for a reactive pose. The box is built
along the actual gorge axis per structure, from the serine Og to the peripheral site
centroid. Gorge geometry differs between species, which is the whole point, so a shared box
would impose one species' geometry on all of them.

*The conversion trap.* Open Babel does not understand HID, and asked to protonate at pH 7.4
it can put the proton on Ne2, undoing the work above silently in every structure. The second
script converts and then inspects the result, checking for a polar hydrogen on Nd1 and its
absence on Ne2. A receptor failing this check is not docked.

---

## 8D dock

```
sbatch --array=0-N stages/8D-dock.slurm RECEPTORS/pdbqt LIGANDS DOCKED
```

**Writes** `DOCKED/<species>/scores.tsv`, `DOCKED/<species>/*.pdbqt`

One array task per receptor, all ligands inside, each receptor using its own box read from
`receptor_prep.tsv`. Exhaustiveness 32 rather than the default 8: the gorge is a deep narrow
channel and the default under-samples it, which shows up as several kJ/mol of run-to-run
scatter on identical input. Three fixed seeds per pair, so docking noise can be quantified
the same way structural noise was.

Only receptors with status `ok` in the tautomer check are docked.

---

## 9A, 9B docking analysis

```
python3 stages/9A-analyse_docking.py --docked DOCKED --receptors RECEPTORS/prepped \
    --prep RECEPTORS/receptor_prep.tsv --aln CLASSIFY/aligned.a2m \
    --endpoints REPORT/analysis_set.tsv --outdir DOCKANALYSIS
python3 stages/9B-analyse_biphasic.py --docked DOCKED \
    --prep RECEPTORS/receptor_prep.tsv --endpoints REPORT/analysis_set.tsv --outdir BIPHASIC
```

**Writes** `DOCKANALYSIS/affinity_species.tsv`, `affinity_variance.tsv`,
`nac_geometry.tsv`, `contacts.tsv`; `BIPHASIC/poses_all.tsv`,
`penetration_descriptors.tsv`, `penetration_icc.tsv`, `biphasic_species.tsv`

Four things, in the order the conclusions depend on each other: affinity reliability across
seeds; near-attack geometry for the reactive ligand, which is physically interpretable in a
way a score is not; contact fingerprints mapped through the alignment so that "touches
position 288" means the same thing in a fish and a water flea; and the control diagnostic,
which correlates each ligand's affinity against the endpoint separately.

The biphasic script re-reads every mode from every seed rather than the single best pose,
projects each onto the gorge axis, and scores the peripheral site, the mid-gorge basin and
the acylation site separately. No new docking is required. The difference between the
peripheral and acylation phases is closer to the quantity that actually varies between
species than any single affinity is.

---

## 10A, 10B model

```
python3 stages/10A-phylo_ssd.py --report REPORT --classify CLASSIFY --tree TREE --outdir FIGURES
python3 stages/10B-model_sensitivity.py --endpoints REPORT/analysis_set.tsv \
    --descriptors STRUCT/descriptors_species.tsv \
    --descriptor-icc STRUCT/variance_partition.tsv \
    --penetration BIPHASIC/penetration_descriptors.tsv \
    --penetration-icc BIPHASIC/penetration_icc.tsv \
    --affinity DOCKANALYSIS/affinity_species.tsv \
    --affinity-icc DOCKANALYSIS/affinity_variance.tsv \
    --geometry DOCKANALYSIS/nac_geometry.tsv \
    --tree <dated species tree> --outdir MODEL
```

**Writes** `FIGURES/fig1..fig6`, `MODEL/univariate_results.tsv`, `MODEL/phylo_signal.tsv`,
`MODEL/variance_partition.tsv`, `MODEL/provenance.tsv`

Run `analyse.py` early, before you invest in docking. Step 4 of it estimates the
phylogenetic signal in the endpoint. If that is near 1, most of the interspecies variance is
already explained by shared ancestry and whatever a structural predictor can add is confined
to the residual. That is a publishable result either way and it costs a minute here against
months of docking.

`model_sensitivity.py` carries the controls that make the rest of the pipeline meaningful:

- **A provenance guard, which is a hard stop.** If a predictor table and its ICC table share
  no descriptor names, they are from different runs and the script refuses to proceed. This
  is not hypothetical; it caught exactly that on the reference dataset.
- **Clade confound decomposition.** For every predictor, eta squared on taxonomic class,
  compared against the endpoint's own. A predictor more taxonomic than the thing it claims
  to explain is the arithmetic signature of a confound, and unlike phylogenetic regression
  this does not depend on a tree.
- **ICC gating.** Disattenuating by one over the square root of a small ICC turns a weak
  correlation into the largest number in the table. Disattenuated values are suppressed
  below `--min-icc`, and predictors with no noise floor at all are marked and excluded
  unless overridden.
- **Sign check.** Under the binding hypothesis, stronger predicted binding should mean
  greater sensitivity. Affinity predictors carry an explicit flag when they do not.
- **Leave-one-out and held-out-clade.** At n around 30 one species is three per cent of the
  data, and withholding a single species from a phylogenetically autocorrelated set leaves
  its relatives in the training partition. Both are reported; the clade-level figure is the
  honest one.

---

## 11 report

```
python3 stages/11-build_report.py --root . --outdir DASHBOARD
```

**Writes** `DASHBOARD/index.html`

One file, data embedded, no server and no network. The validation status is computed before
anything is displayed and stated in the banner. If no predictor survived correction, it says
so at the top.

What it projects for unmeasured species is not the structural descriptors. When phylogenetic
signal in the endpoint is high, the best available predictor of an untested species is its
relatives, so prediction is phylogenetic imputation under Brownian motion with credible
intervals that widen with distance from the nearest measured relative. That is read-across,
made quantitative, and it is honest about what the data supports.

---

## Scripts not yet released

Four scripts are referenced by the driver and are not in this first release. Two of them
are numbered stages (2 and 4A) and appear in `bin/sentinel stages` marked NOT IN THIS
RELEASE. `bin/sentinel`
stops with an explanation when it reaches one.

**`2-retrieve_orthologues.py`** (stage 2) queries OrthoDB for the orthologous group, writes one
FASTA per species into `ODB/by_species/`, and builds `species_inventory.tsv` with the NCBI
taxid and taxonomic class per species plus `ecotox_coverage.tsv` mapping endpoint species
onto available sequences, including genus-level proxies. Substituting it means producing
those three things by any means you like; nothing downstream cares how they were made.

**`4A-finalise_endpoints.py`** (stage 4A) collapses the cleaned ECOTOX table to one row per species,
resolves congeners onto a shared sequence proxy by geometric mean, and writes
`FINAL/endpoint_species.tsv` with the OrthoDB tags each endpoint species maps to.

**`make_dataset_report.py`** writes `REPORT/analysis_set.tsv`, which is the endpoint table
used by stages 9A and 10B. In its absence the driver falls back to
`FINAL_SET/endpoint_final.tsv`, which carries the same species and endpoints under different
column names, so check the column mapping before trusting the fallback.

**`timetree_helper.py`** reconciles the species set against a dated TimeTree export and
writes the pruned species tree used for the phylogenetic comparative models. Without it, set
`SPECIES_TREE` to a Newick file you have pruned yourself, or accept the gene-tree fallback
and its caveats.

Also absent: `setup_dock_env.sh`, which built the micromamba environment holding vina,
obabel, Meeko and rdkit on the reference system, and the figure scripts
(`make_figures.py`, `plot_composite.py`) used for the manuscript panels.
