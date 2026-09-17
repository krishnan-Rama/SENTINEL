# SENTINEL

**Species Sensitivity from Target Interaction and Evolutionary Lineage**

SENTINEL models cross-species chemical susceptibility by tracing evolutionary divergence in xenobiotic target proteins. By integrating comparative sequence analysis, 3D structural modelling, and binding simulations, the workflow translates molecular initiating events (MIEs) into predictive Species Sensitivity Distributions (SSDs) across uncharacterised taxa.

The pipeline is engineered for SLURM-managed HPC clusters, distributing GPU-accelerated structural predictions and parallel CPU tasks across cluster nodes.

---

## Pipeline Architecture

```mermaid
flowchart TD
    subgraph INP["<h3>Input Data<h3>"]
        direction TB
        I1["Target Sequences / Accessions"]
        I2["Chemical Structures (.sdf / SMILES)"]
        I3["Empirical Ecotoxicity Endpoints (LC50/EC50)"]
    end

    subgraph ST1["<h3>Phase 1: Orthology & Lineage [SLURM CPU]</h3>"]
        direction TB
        S1A["Target Homology & Domain Slicing<br/><b>HMMER v3.4</b> • <b>MAFFT v7.5</b>"]
        S1B{"Tree Engine"}
        S1C["Lineage Phylogeny<br/><b>FastTree</b> (Default)"]
        S1D["Bootstrapped ML Tree<br/><b>RAxML-NG</b> (--with-ml-tree)"]

        S1A --> S1B
        S1B -->|Default| S1C
        S1B -->|Optional| S1D
    end

    subgraph ST2["<h3>Phase 2: Structural Modelling [SLURM GPU]<h3>"]
        direction TB
        S2A["Ensemble Folding & Quality Filtering<br/><b>LocalColabFold (AlphaFold2)</b>"]
        S2B["Explicit-Solvent Minimisation<br/><b>GROMACS</b> (--with-minimise)"]

        S2A -.->|Optional| S2B
    end

    subgraph ST3["<h3>Phase 3: Molecular Docking [SLURM CPU]<h3>"]
        direction TB
        S3A["Ligand & Receptor Preparation<br/><b>RDKit</b> • <b>Meeko</b> • <b>Open Babel</b>"]
        S3B["Conformational Docking Arrays<br/><b>AutoDock Vina</b>"]

        S3A --> S3B
    end

    subgraph ST4["<h3>Phase 4: Evolutionary Sensitivity [Master Node]<h3>"]
        direction TB
        S4A["Evolutionary Regressions & Sensitivity<br/><b>PGLS</b> • <b>Scikit-learn</b>"]
        S4B["Phylogenetic Signal & Validation<br/><b>Pagel's &lambda;</b> • <b>Clade Cross-Validation</b>"]

        S4A --> S4B
    end

    subgraph OUT["<h3>Primary Outputs<h3>"]
        direction TB
        O1["master_table.csv<br/><i>(Alignments, calls & metadata)</i>"]
        O2["MODEL/univariate_results.tsv<br/><i>(Signal & clade error metrics)</i>"]
        O3["DASHBOARD/index.html<br/><i>(Interactive report & SSD curves)</i>"]
    end

    %% Wiring Inputs to Pipeline
    I1 -->|--entry sequences| S1A
    S1C -->|--entry structures| S2A
    S1D -->|--entry structures| S2A
    S2A -->|--entry poses| S3A
    S2B -->|--entry poses| S3A
    I2 --> S3A
    S3B -->|--entry analysis| S4A
    I3 --> S4A

    %% Outputs
    S4B --> O1
    S4B --> O2
    S4B --> O3

    %% High-Contrast Styling (WCAG AAA Compliant)
    classDef io fill:#1e293b,stroke:#0f172a,stroke-width:2px,color:#ffffff;
    classDef cpu fill:#14532d,stroke:#052e16,stroke-width:2px,color:#ffffff;
    classDef gpu fill:#7f1d1d,stroke:#450a0a,stroke-width:2px,color:#ffffff;
    classDef master fill:#1e3a8a,stroke:#172554,stroke-width:2px,color:#ffffff;

    class I1,I2,I3,O1,O2,O3 io;
    class S1A,S1B,S1C,S1D,S3A,S3B cpu;
    class S2A,S2B gpu;
    class S4A,S4B master;
```

---

## Quickstart

### 1. Environment Setup

```bash
git clone https://github.com/krishnan-Rama/SENTINEL.git
cd SENTINEL

conda env create -f environment.yml
conda activate sentinel
```

Verify binary discovery, cluster modules, and CLI accessibility:

```bash
bin/sentinel --help
bin/sentinel doctor
```

### 2. Configuration and Submission

All target specifications, active-site boundaries, ligand paths, and compute parameters are configured in a project environment file.

```bash
# Prepare a run profile from the template
cp config/project.env config/my_run.env
$EDITOR config/my_run.env

# Validate inputs and inspect the execution plan
bin/sentinel -c config/my_run.env -n run

# Dispatch array jobs to SLURM
bin/sentinel -c config/my_run.env \
    --hpc slurm \
    --partition standard \
    --gpu-partition gpu \
    --account proj123 \
    run
```

---

## Command-Line Interface

```text
bin/sentinel -c <config.env> [OPTIONS] <COMMAND>
```

### Commands

* `run`: Execute the planned pipeline workflow.
* `doctor`: Audit environment modules, external binaries, and input file integrity.
* `--help`: Display the CLI syntax and default arguments.

### Pipeline Execution Options

* `-n`: Dry-run mode. Displays planned execution steps and task counts without dispatching jobs.
* `--entry <STAGE>`: Resume or start from intermediate checkpoints.
* `--force`: Disregard cached outputs and force recomputation of downstream stages.
* `--with-ml-tree`: Infer a maximum-likelihood phylogeny via RAxML-NG instead of FastTree.
* `--with-minimise`: Execute explicit-solvent energy minimisation via GROMACS prior to docking.

### SLURM Resource Flags

* `--hpc slurm|local`: Execution mode (default: `local`).
* `--partition <NAME>`: CPU partition for sequence alignment, tree inference, and docking.
* `--gpu-partition <NAME>`: GPU partition for ColabFold structural predictions.
* `--account <NAME>`: SLURM project allocation code.
* `--max-parallel <N>`: Maximum number of concurrent array tasks.
* `--threads <N>`: Number of CPU threads per task.

---

## Modular Entry Points

Use `--entry` to resume from intermediate checkpoints or inject existing structural data. Completed stages are cached automatically.

| Current Data State | Command | Scope of Execution |
|---|---|---|
| Species accessions and endpoints | `--entry sequences` (default) | Executes full pipeline from alignment onward. |
| Precomputed 3D structures | `--entry structures --models <DIR> --alignment <FILE>` | Skips ColabFold; proceeds to docking. |
| Docked receptor-ligand poses | `--entry poses --docked <DIR> --alignment <FILE>` | Skips docking; extracts descriptors and fits models. |
| Feature matrices and affinities | `--entry analysis` | Executes regressions, phylogenetic models, and report generation. |

---

## Primary Outputs

All generated artefacts are written to the directory specified by `--outdir` in your configuration file:

* `master_table.csv`: Consolidated dataset of sequence accessions, active-site residue calls, orthologue validation, and harmonised endpoints.
* `MODEL/univariate_results.tsv`: Statistical summaries per descriptor, including phylogenetic signal (Pagel's lambda) and leave-one-clade-out cross-validation error.
* `DASHBOARD/index.html`: Self-contained interactive report displaying structural alignments, binding poses, and predicted SSD curves.

---

## Contact

Rama Krishnan  
School of Biosciences, Cardiff University  
Email: krishnanr1@cardiff.ac.uk
