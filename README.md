# SENTINEL

**S**pecies s**EN**sitivity from **T**arget **IN**teraction and **E**volutionary **L**ineage.

SENTINEL takes a molecular initiating event, a chemical, and a list of species, and returns
an estimate of how sensitivity to that chemical is distributed across those species,
together with the evidence for whether the estimate should be believed.

It is built for regulatory scientists and ecotoxicologists who need to say something
defensible about a species that has never been tested. It is not a docking front-end and
not a structure-prediction wrapper, although it uses both. The part that matters is the
audit trail: at every stage it records why a species was kept or dropped, how much of a
measurement is signal and how much is noise, and whether a correlation survives the
controls. Where the evidence does not support a projection, it says so instead of returning
a number.

The reference application is acetylcholinesterase and the organophosphate insecticide
chlorpyrifos across 81 metazoan species, which is the dataset the accompanying paper
reports.

> **Status: working research code.** It runs, it has produced a complete dataset, and the
> outputs are what the paper is based on. It is not yet portable to an arbitrary target
> without editing Python, and four helper scripts are not in this first release. See
> [Limitations](#limitations).

---

## Install and run

```bash
git clone https://github.com/<you>/sentinel.git
cd sentinel

conda env create -f environment.yml && conda activate sentinel
# or: pip install -r requirements.txt, plus the external tools listed below

cp config/project.env config/my_run.env
$EDITOR config/my_run.env

bin/sentinel -c config/my_run.env doctor     # check tools, inputs, what is missing
bin/sentinel -c config/my_run.env -n run     # print the plan without running it
bin/sentinel -c config/my_run.env run        # go
```

`bin/sentinel --help` is the full reference. Four commands: `run`, `status`, `stages`,
`doctor`.

---

## Starting where you already are

Most people arrive with part of the work already done. `--entry` sets the starting point,
and the driver checks that you have supplied what that entry point needs.

| You have | Command |
|----------|---------|
| A species list and a toxicity table | `--entry sequences` (the default) |
| Predicted structures, no docking | `--entry structures --models DIR --alignment FILE` |
| Docking results already | `--entry poses --docked DIR --alignment FILE` |
| Descriptor and affinity tables | `--entry analysis` |

```bash
# full run on a cluster
bin/sentinel -o ~/runs/cpf --gpu-partition gpu --partition epyc --account proj123 run

# AlphaFold models in hand, take it from quality control onwards
bin/sentinel -o ~/runs/cpf --entry structures \
    --models ~/af_models --alignment ~/runs/cpf/CLASSIFY/aligned.a2m run

# orthologue verification and structures only, no docking at all
bin/sentinel -o ~/runs/cpf --no-dock run

# redo just the association testing with a proper dated tree
bin/sentinel -o ~/runs/cpf --species-tree timetree.nwk --only 10B,11 run
```

Finer control than the presets: `--from 5A`, `--to 8D`, `--only 8A,8B,8C`, `--skip 4C`.
Optional stages are off unless asked for: `--with-ml-tree` adds the bootstrapped publication
phylogeny, `--with-minimise` adds explicit-solvent minimisation.

A stage whose key output already exists is skipped and reported as cached, so an interrupted
run resumes by rerunning the same command. `--force` overrides that.

One caveat on `--entry structures`. You still need a profile alignment, because the residue
numbering frame is what tells the pipeline which residue in your model is the catalytic one.
Either point `--alignment` at an existing A2M or run stages 1 to 3A first. The driver refuses
rather than guessing.

---

## Cluster options

Everything scheduler-related is a flag, so nothing inside the SLURM scripts needs editing
for a different site:

```
--hpc slurm|local      submit array jobs, or just print the sbatch commands
--partition NAME       CPU partition: docking, phylogeny, minimisation
--gpu-partition NAME   GPU partition: structure prediction
--account NAME         --qos NAME
--time-fold / --time-dock / --time-tree
--max-parallel N       cap concurrent array tasks
--threads N            --chunk N
```

`--hpc local` prints the exact `sbatch` line for each heavy stage instead of submitting it,
which is also what you want if your site uses PBS or SGE and you need to translate.

Array sizes are computed, not typed. Stage 5B sizes itself from the number of folding
chunks; stage 8D sizes itself from the number of receptors that passed the tautomer check.

---

## The stages

```
1     verify the reference backbone by accession
2     retrieve candidate orthologues per species
3A    positional active-site classification
3B    clade assignment by tree, one representative per species
3C    bootstrapped publication phylogeny                  (optional)
4A    collapse endpoints to one row per species
4B    build the endpoint and prediction sets
4C    master table, one row per sequence
5A    trim to the catalytic domain, write folding chunks
5B    predict five structural models per sequence         (GPU array job)
5C    explicit-solvent minimisation                       (optional)
6     active-site pLDDT and the between-seed noise floor
7     gorge descriptors on all models, with ICC
8A    build and identity-check the ligand panel
8B    catalytic geometry and a per-species docking box
8C    PDBQT conversion with a tautomer check
8D    dock, three seeds per receptor-ligand pair          (CPU array job)
9A    affinity reliability, near-attack geometry, contacts
9B    biphasic site occupancy along the gorge axis
10A   phylogenetic signal, SSD, per-step figures
10B   association testing, confound decomposition, validation
11    single-file gated HTML report
```

Scripts live in `stages/` under exactly these names, so the numbering is the order and the
filename is the stage. Each is a standalone script with a documented argument list; the
driver only resolves paths, checks inputs and dispatches jobs. `-n` prints the command it
would run so you can copy it.

Only 5B and 8D are expensive, and both are parallel across species. On the reference system
the 81-species dataset took roughly 30 hours of wall clock for folding across seven
concurrent array tasks and roughly 20 hours for docking across 82. Both depend heavily on
queue behaviour. Treat them as an order of magnitude and record your own.

Per-stage inputs, outputs and rationale: [docs/STAGES.md](docs/STAGES.md).

---

## What it actually does

Six things, in order, each of which can fail loudly.

**Establishes what the target is.** Reference sequences are fetched by UniProt accession
rather than by gene-name search, and each is checked against an expected length window and a
required substring of the protein name. This is not defensive programming for its own sake:
in an earlier version a gene-name query for human carboxylesterase 1 silently returned
metallothionein-2, a 61-residue protein, which would have rooted the whole phylogeny on
nonsense without any error being raised.

**Finds the orthologue in every species, and proves it.** Candidates are aligned to a profile
built from the verified references, scored at diagnostic active-site positions, then placed
on a maximum-likelihood tree. Clade membership decides what a sequence is; the residues are
reported as evidence. In the reference dataset this corrected 56 sequences that residue
heuristics alone had called acetylcholinesterase but which are butyrylcholinesterases,
carboxylesterases or catalytically dead family members.

**Audits the applicability domain.** Attrition is tracked from the source toxicity database
to the final analysis set, and the pipeline tests whether the species that were lost differ
in sensitivity from the species that were kept. A screening tool whose training set is
biased toward tolerant taxa will be confidently wrong, and there is no way to detect that
after the fact.

**Measures its own noise floor.** Every target is predicted five times from different seeds
and every descriptor is computed on all five. Variance is partitioned into between-species
and between-seed components, and descriptors whose seed scatter rivals their species range
are excluded before any statistics are run. In the reference dataset four of eleven
descriptors failed this check, including the one with the most intuitive mechanistic story
attached to it.

**Docks negative controls on purpose.** Alongside the active compound the panel includes
molecules chemically incapable of the reaction the mechanism depends on. If an inert
molecule ranks species as well as the active one does, the ranking is measuring generic
pocket accommodation, and the pipeline reports that rather than hiding it.

**Separates target signal from ancestry.** Species are not independent observations. For
every predictor the pipeline reports how much of its variance sits between taxonomic
classes, refits the correlation inside each clade, runs phylogenetic regression, and
cross-validates by withholding whole clades rather than single species. On the reference
dataset this is what turned an apparently strong structure-toxicity relationship into what
it actually is: a phylogenetic prior on baseline sensitivity.

---

## What you get

Everything is written under `--outdir` as tab-separated tables, so any stage can be
inspected or replaced without re-running the ones before it.

`master_table.csv` has one row per retrieved sequence, not one per species, and keeps the
butyrylcholinesterases, carboxylesterases and non-catalytic family members that the target
calls were discriminated against. A reviewer who wants to check a classification can do so
without re-running anything.

`MODEL/univariate_results.tsv` has one row per predictor with its correlation against the
endpoint, its reliability, how much of its variance is taxonomic, its correlation refitted
inside each clade, its phylogenetic regression and its held-out-clade error. This is the
table the paper's conclusions rest on.

`DASHBOARD/index.html` is a single file with the data embedded, no server and no network
access, so it can be emailed or archived alongside a dossier. It is gated: validation status
is computed first and stated in the banner, so a ranking is never displayed without the
evidence for whether it means anything.

---

## External tools

| Tool | Stage | Notes |
|------|-------|-------|
| MAFFT 7.5 or later | 3A, 3B, 3C | L-INS-i for the publication tree; `--auto` silently degrades at this scale |
| HMMER 3.4 or later | 3A | profile alignment |
| FastTree | 3B | fast tree for clade assignment |
| RAxML-NG 1.2 or later | 3C | publication tree with bootstraps |
| LocalColabFold | 5B | needs a GPU, and `hhsearch` on PATH if templates are enabled |
| AutoDock Vina | 8D | |
| Open Babel, Meeko | 8A, 8C | PDBQT conversion |
| RDKit | 8A | ligand construction and identity checks |
| GROMACS | 5C | only if you go on to molecular dynamics |

`bin/sentinel doctor` tells you which of these it can find.

---

## Applying it to a different target

SENTINEL is target-agnostic in design and target-specific in this release. Porting it needs
four things, three of which are edits to a Python file rather than a config change:

A **verified reference set** for the new protein family, including catalytically inactive
outgroups. Edit `ANCHORS` in `stages/1-anchors.py`. Every entry needs an accession, a length
window and a name substring, and the script refuses to proceed if any fails.

A **residue numbering frame** from an experimental structure, and a motif that locates it.
The current code finds the catalytic serine through the nucleophile elbow motif and anchors
everything to Torpedo mature numbering. Edit `TORPEDO_POS` and the motif regular expression
in `stages/3A-classify_active_site.py`, then set `--domain START END` for the new domain
boundaries.

A **ligand panel with inactive analogues.** Edit `PANEL` in `stages/8A-prepare_ligands.py`.
The controls are the point: without a molecule that cannot perform the mechanism, the
specificity check does not exist.

A **cleaned toxicity table.** Format in [docs/INPUTS.md](docs/INPUTS.md).

Making the first three configuration rather than code is the main item for the next version.

---

## Limitations

Read this before relying on anything it produces.

**Not all of the workflow is in this release.** Stages 2 and 4A are not included, along with
the dataset report and the TimeTree helper. `bin/sentinel` stops with an explanation rather
than a traceback when it reaches one, and `doctor` lists them. What each does, and what you
can substitute, is in [docs/STAGES.md](docs/STAGES.md#scripts-not-yet-released).

**The default tree is a gene tree, and that is not good enough for the final analysis.** A
gene tree of a duplicated family groups by paralogue before taxonomy and its branch lengths
are divergence rather than time. It is fine for deciding what a sequence is, which is stage
3B. It is not fine for phylogenetic regression or read-across. Pass `--species-tree`. If you
do not, stage 10B warns and continues. Two runs of the reference dataset differing only in
this gave Pagel's lambda of 0.885 and 0.965, which is enough to change conclusions.

**Docking scores are not potencies.** For a covalent inhibitor, potency depends on both the
binding of the pre-reaction complex and the barrier to the reaction itself, and a docking
score touches only the first. On the reference dataset no affinity predicted toxicity, and
five of six ran in the direction opposite to the mechanistic expectation. The informative
comparisons are between ligands, not the absolute numbers.

**Rigid receptors.** One relaxed model per species, no receptor flexibility. Search
stochasticity is quantified across three seeds; conformational uncertainty is not.

**Small n.** Structural extrapolation is only testable against species that have both a
target and a measured endpoint. In the reference dataset that was 30 species out of 84
starting points, which supports detection of moderate to large effects and nothing subtler.

---

## Citing

Please cite the paper (details on acceptance) and this repository. `CITATION.cff` is
included and GitHub renders it as a citation box.

## Licence

MIT.

## Contact

Rama Krishnan, School of Biosciences, Cardiff University. Issues and pull requests welcome.
If a stage failed, the output of `bin/sentinel doctor` plus the relevant `*_report.tsv` or
`*_qc.tsv` is usually enough to diagnose it.
