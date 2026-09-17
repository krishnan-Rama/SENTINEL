# SENTINEL

**S**pecies s**EN**sitivity from **T**arget **IN**teraction and **E**volutionary **L**ineage

SENTINEL automates the derivation and prediction of cross-species Species Sensitivity Distributions (SSDs). It bridges molecular initiating events and species-level toxicity thresholds by coupling evolutionary phylogenetics, AlphaFold structural ensembles, and molecular docking affinities into an auditable modelling pipeline.

---

## Requirements

* **OS:** Linux (x86_64)
* **Environment:** Conda or Mamba
* **Hardware:** CUDA-capable GPU (required for ColabFold structural prediction; not required if using `--entry structures` or precomputed models)
* **Cluster:** SLURM workload manager (optional; supported natively for parallel array jobs)

---

## Installation

```bash
# Clone the repository
git clone https://github.com/krishnan-Rama/SENTINEL.git
cd SENTINEL

# Set up Conda environment
conda env create -f environment.yml
conda activate sentinel

# Verify environment, paths, and external binaries
bin/sentinel doctor
```

> **Note:** External binaries such as `LocalColabFold` and `GROMACS` must reside in your system `$PATH` if running structural inference or explicit-solvent relaxation.

---

## Quickstart

### 1. Set Up Configuration

Copy the template configuration file and define your run parameters:

```bash
cp config/project.env config/my_run.env
$EDITOR config/my_run.env
```

Specify your target protein, reference UniProt anchors, catalytic boundaries, and ligand structures directly inside `config/my_run.env`.

### 2. Validate and Plan

```bash
# Audit inputs, tools, and path availability
bin/sentinel -c config/my_run.env doctor

# Dry-run execution plan without launching jobs
bin/sentinel -c config/my_run.env -n run
```

### 3. Run Pipeline

Execute locally:
```bash
bin/sentinel -c config/my_run.env run
```

Or dispatch heavy array jobs to a SLURM cluster:
```bash
bin/sentinel -c config/my_run.env \
    --hpc slurm \
    --partition standard \
    --gpu-partition gpu \
    --account proj123 \
    run
```

---

## Modular Checkpoints & Reusable Metadata

SENTINEL supports resuming from intermediate stages or executing downstream statistical models directly using provided metadata tables:

| Input Available | Entry Command | Description |
|---|---|---|
| Species accessions & raw data | `--entry sequences` | Full pipeline execution (Default). |
| Precomputed 3D structures | `--entry structures --models DIR --alignment FILE` | Bypasses sequence search, MSA, and ColabFold. |
| Docked receptor-ligand complexes | `--entry poses --docked DIR --alignment FILE` | Bypasses docking; extracts interaction descriptors. |
| Curated metadata & affinity tables | `--entry analysis` | Skips all structural work; runs phylogenetic regressions and SSDs. |

### Standalone Metadata Assets

If you do not need to recompute structures or docking poses, you can use or adapt the processed metadata files directly:

* **`endpoint_final.tsv`**: Standardised and harmonised species toxicity endpoints across tested taxa. Useful for training independent machine learning models, re-fitting SSD curves, or validating species-specific thresholds.
* **`prediction_final.tsv`**: Compiled docking descriptors, binding affinities, active-site classifications, and predicted species sensitivities. Adaptable for custom statistical pipelines, phylogenetic generalized least squares (PGLS), or exploratory data analysis in R/Python.

To execute evolutionary regressions and dashboard generation directly from these tables:
```bash
bin/sentinel -c config/my_run.env --entry analysis --endpoints endpoint_final.tsv --predictions prediction_final.tsv run
```

---

## Execution Options

* `--force`: Force recalculation of stages, bypassing cached intermediates.
* `--with-ml-tree`: Infer maximum-likelihood phylogeny via RAxML-NG (replaces FastTree default).
* `--with-minimise`: Run explicit-solvent energy minimisation via GROMACS prior to docking.
* `--max-parallel N`: Set concurrency limit for cluster array jobs.

---

## Pipeline Workflow

```
Raw Sequences & Chemical Structures
              │
              ▼
[1] Orthology & Active-Site Screening  (MAFFT / HMMER)
              │
              ▼
[2] 3D Structural Ensembles            (LocalColabFold)
              │
              ▼
[3] Receptor Preparation & Docking     (AutoDock Vina / Meeko)
              │
              ▼
[4] Evolutionary Sensitivity Models    (PGLS / Machine Learning)
              │
              ▼
[5] Interactive Validation Dashboard   (HTML / Standalone)
```

---

## Primary Outputs

All generated artefacts are written to your configured output directory (`--outdir`):

* `master_table.csv`: Curated sequence metadata, active-site residue calls, and aligned annotations.
* `MODEL/univariate_results.tsv`: Model fit summaries, Pagel's $\lambda$ phylogenetic signals, and leave-one-clade-out cross-validation error metrics.
* `DASHBOARD/index.html`: Self-contained interactive report featuring structural superpositions, predicted SSD curves, and diagnostic plots for offline sharing.

---

## Input Formats

For detailed column schemas, sequence header formats, and chemical structure requirements (SDF/MOL2), see [docs/INPUTS.md](docs/INPUTS.md).

---

## Contact

**Rama Krishnan**  
School of Biosciences, Cardiff University  
Email: krishnanr1@cardiff.ac.uk
