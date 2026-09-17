# SENTINEL

**Species Sensitivity from Target Interaction and Evolutionary Lineage**

SENTINEL models cross-species chemical susceptibility by tracing evolutionary divergence in xenobiotic target proteins. By integrating comparative sequence analysis, 3D structural modelling, and binding simulations, the workflow translates molecular initiating events (MIEs) into predictive Species Sensitivity Distributions (SSDs) across uncharacterised taxa.

The pipeline is engineered for SLURM-managed HPC clusters, distributing GPU-accelerated structural predictions and parallel CPU tasks across cluster nodes.

---

## Quickstart

### 1. Environment Setup

```bash
git clone https://github.com/krishnan-Rama/SENTINEL.git
cd SENTINEL

conda env create -f environment.yml
conda activate sentinel
```

### 2. Configuration and Submission

All run parameters, sequence identifiers, domain boundaries, and ligand structures are configured in a run environment file.

```bash
# Generate and edit the run profile for your target system
cp config/project.env config/my_run.env
$EDITOR config/my_run.env

# Dry run: validate inputs and display the execution plan
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

## Command-Line Interface

Because the executable does not provide an interactive `--help` menu, use the invocation syntax and options documented below:

```text
bin/sentinel -c <config.env> [OPTIONS] run
```

### Execution Options

* `-n`: Dry-run mode. Displays planned execution steps and task counts without dispatching jobs.
* `--entry <STAGE>`: Resume or start from intermediate checkpoints.
  * `sequences`: Default. Runs the full workflow from initial sequence sets.
  * `structures --models <DIR> --alignment <FILE>`: Skips structure prediction; proceeds to docking.
  * `poses --docked <DIR> --alignment <FILE>`: Skips docking; extracts descriptors and fits models.
  * `analysis`: Runs only phylogenetic models, cross-validation, and report generation.
* `--force`: Ignores cached outputs and forces recomputation of specified stages.
* `--with-ml-tree`: Infers a maximum-likelihood phylogeny using RAxML-NG instead of the FastTree default.
* `--with-minimise`: Runs explicit-solvent energy minimisation via GROMACS prior to docking.

### SLURM Flags

* `--hpc slurm|local`: Execution mode (default: `local`).
* `--partition <NAME>`: CPU partition for alignments, tree inference, and docking.
* `--gpu-partition <NAME>`: GPU partition for ColabFold structural predictions.
* `--account <NAME>`: SLURM allocation or project account code.
* `--max-parallel <N>`: Maximum number of concurrent array tasks.
* `--threads <N>`: Number of CPU threads per task.

---

## Primary Outputs

All pipeline artefacts are written to the directory specified by `--outdir` in your configuration file:

* `master_table.csv`: Consolidated dataset of sequence accessions, active-site residue calls, orthologue validation, and harmonised endpoints.
* `MODEL/univariate_results.tsv`: Statistical summaries per descriptor, including phylogenetic signal (Pagel's $\lambda$) and leave-one-clade-out cross-validation error.
* `DASHBOARD/index.html`: Self-contained interactive report displaying alignments, binding poses, and predicted SSD curves.

---

## Contact

Rama Krishnan  
School of Biosciences, Cardiff University  
Email: `krishnanr1@cardiff.ac.uk`
