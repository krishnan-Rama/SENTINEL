# SENTINEL

**S**pecies s**EN**sitivity from **T**arget **IN**teraction and **E**volutionary **L**ineage

SENTINEL predicts cross-species chemical sensitivity distributions (SSDs) using target-ligand interactions and evolutionary phylogenetics. It provides an automated, auditable pipeline linking molecular initiating events to species-level sensitivity endpoints.

---

## Installation

```bash
git clone https://github.com/krishnan-Rama/SENTINEL.git
cd sentinel

conda env create -f environment.yml && conda activate sentinel
```

Run `bin/sentinel --help` for the complete command-line reference.

---

## Quickstart

```bash
# 1. Configure run parameters
cp config/project.env config/my_run.env
$EDITOR config/my_run.env

# 2. Check environment and inputs
bin/sentinel -c config/my_run.env doctor

# 3. Dry run: display execution plan
bin/sentinel -c config/my_run.env -n run

# 4. Execute pipeline
bin/sentinel -c config/my_run.env run
```

---

## Pipeline Entry Points

Resume or start the pipeline from intermediate data using `--entry`:

| Existing Data | Command |
|---|---|
| Species list and toxicity table | `--entry sequences` (default) |
| Predicted structures | `--entry structures --models DIR --alignment FILE` |
| Docking poses | `--entry poses --docked DIR --alignment FILE` |
| Descriptors and affinities | `--entry analysis` |

### Execution Controls

* **Stage selection:** `--from 5A`, `--to 8D`, `--only 8A,8B,8C`, or `--skip 4C`.
* **Optional stages:** `--with-ml-tree` (bootstrapped publication phylogeny) and `--with-minimise` (explicit-solvent energy minimisation).
* **Caching:** Completed stages are automatically skipped on rerun. Use `--force` to recompute.

```bash
# Example: SLURM submission starting from AlphaFold models
bin/sentinel -o ~/runs/cpf --entry structures \
    --models ~/af_models --alignment ~/runs/cpf/CLASSIFY/aligned.a2m \
    --gpu-partition gpu --partition epyc --account proj123 run
```

---

## HPC and Cluster Configuration

Job arrays for structural prediction (stage 5B) and docking (stage 8D) scale automatically to the input size.

| Flag | Description |
|---|---|
| `--hpc slurm\|local` | Submit array jobs via SLURM or output shell commands locally |
| `--partition NAME` | CPU partition for docking, tree inference, and minimisation |
| `--gpu-partition NAME` | GPU partition for ColabFold structure prediction |
| `--account NAME` | SLURM account or allocation identifier |
| `--max-parallel N` | Maximum number of concurrent array tasks |
| `--threads N` | Thread allocation per task |

---

## Pipeline Stages

The workflow organises modular scripts into six core operational phases:

| Phase | Stages | Description |
|---|---|---|
| **1. Orthology & Phylogeny** | `1` to `3C` | Validate reference sequences, classify active-site residues, and infer lineage trees. |
| **2. Endpoint Assembly** | `4A` to `4C` | Harmonise toxicity endpoints across species and compile master dataset tables. |
| **3. Structural Modelling** | `5A` to `7` | Predict 3D model ensembles, screen pLDDT variance, and extract pocket descriptors. |
| **4. Molecular Docking** | `8A` to `9B` | Prepare receptor grids, dock target and control ligands, and profile binding poses. |
| **5. Evolutionary Analysis** | `10A` to `10B` | Quantify phylogenetic signal, fit clade-level regressions, and cross-validate. |
| **6. Reporting** | `11` | Compile metrics into an interactive, self-contained HTML dashboard. |

Further details on per-stage inputs and outputs are documented in [docs/STAGES.md](docs/STAGES.md).

---

## Outputs

All artefacts are saved to `--outdir`:

* `master_table.csv`: Comprehensive sequence metadata, active-site residue calls, and orthologue classifications.
* `MODEL/univariate_results.tsv`: Statistical summaries per descriptor, including phylogenetic signal, clade-specific correlations, and held-out clade cross-validation errors.
* `DASHBOARD/index.html`: Self-contained interactive report with embedded results and validation summaries for offline sharing and archiving.

---

## External Dependencies

Verify installation paths using `bin/sentinel doctor`:

| Tool | Pipeline Stage | Role |
|---|---|---|
| MAFFT (>= 7.5) | 3A, 3B, 3C | Multiple sequence alignment |
| HMMER (>= 3.4) | 3A | Profile alignment |
| FastTree | 3B | Clade assignment |
| RAxML-NG (>= 1.2) | 3C | Maximum-likelihood phylogenetic inference |
| LocalColabFold | 5B | Structural modelling |
| AutoDock Vina | 8D | Molecular docking |
| Open Babel & Meeko | 8A, 8C | Ligand preparation and PDBQT conversion |
| RDKit | 8A | Chemoinformatics and validation |
| GROMACS | 5C | Energy minimisation (optional) |

---

## Adapting to New Targets

1. **Define target anchors:** Update `ANCHORS` in `stages/1-anchors.py` with UniProt accessions, expected sequence lengths, and naming criteria.
2. **Set catalytic framework:** Update the catalytic motif pattern and reference numbering in `stages/3A-classify_active_site.py`, and supply domain boundaries using `--domain START END`.
3. **Configure ligand panel:** Specify active compounds and inactive chemical controls in `stages/8A-prepare_ligands.py`.
4. **Supply toxicity data:** Format your experimental endpoint table as described in [docs/INPUTS.md](docs/INPUTS.md).

---

## Contact

Rama Krishnan, School of Biosciences, Cardiff University: krishnanr1@cardiff.ac.uk
