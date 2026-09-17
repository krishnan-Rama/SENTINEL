# SENTINEL

**Species Sensitivity from Target Interaction and Evolutionary Lineage**

SENTINEL models cross-species chemical susceptibilities by tracing evolutionary divergence in xenobiotic target proteins. By integrating comparative sequence analysis, 3D structural modelling, and binding simulations, the workflow translates molecular initiating events (MIE) into predictive Species Sensitivity Distributions (SSDs) across uncharacterised taxa.

The pipeline is engineered for SLURM-managed HPC clusters, distributing GPU-accelerated structural predictions and parallel CPU tasks across cluster nodes.

---

## Quickstart

### 1. Environment Setup

```bash
git clone https://github.com/krishnan-Rama/SENTINEL.git
cd SENTINEL

conda env create -f environment.yml
conda activate sentinel

# Audit cluster paths, modules, and external binaries
```

### 2. Configuration and Submission

```bash
# Generate and edit the run profile for your target system
cp config/project.env config/my_run.env
$EDITOR config/my_run.env

# Validate input integrity and inspect the execution plan
bin/sentinel -c config/my_run.env -n run

# Submit array jobs to SLURM
bin/sentinel -c config/my_run.env \
    --hpc slurm \
    --partition standard \
    --gpu-partition gpu \
    --account proj123 \
    run
```

---

## Modular Execution Checkpoints

Resume or start from intermediate stages using `--entry`. Completed steps are cached automatically.

| Current Data State | Command | Scope of Execution |
|---|---|---|
| Species accessions + endpoints | `--entry sequences` (default) | Executes full pipeline from alignment onward. |
| Precomputed 3D structures | `--entry structures --models DIR --alignment FILE` | Skips ColabFold; proceeds to grid prep and docking. |
| Docked receptor-ligand poses | `--entry poses --docked DIR --alignment FILE` | Skips docking; extracts descriptors and runs models. |
| Feature matrices and affinities | `--entry analysis` | Executes regressions, phylogenetic models, and report. |

### Core Flags

* `--force`: Disregards existing cached artefacts and recomputes downstream stages.
* `--with-ml-tree`: Infers a maximum-likelihood phylogeny using RAxML-NG (replaces default FastTree).
* `--with-minimise`: Executes explicit-solvent energy minimisation via GROMACS prior to docking.
* `--max-parallel N`: Caps the number of concurrent SLURM array tasks.

---

## Primary Outputs

All generated artefacts are written to the directory assigned in `--outdir`:

* `master_table.csv`: Consolidated matrix of sequence identifiers, active-site calls, orthology tiering, and harmonised endpoints.
* `MODEL/univariate_results.tsv`: Statistical summaries, phylogenetic signal scores (Pagel's $\lambda$), and leave-one-clade-out cross-validation metrics.
* `DASHBOARD/index.html`: Standalone interactive HTML report compiling structural alignments, docking conformations, and predicted SSD curves.

---

## Contact

Rama Krishnan  
School of Biosciences, Cardiff University  
Email: `krishnanr1@cardiff.ac.uk`
