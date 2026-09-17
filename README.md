# SENTINEL

**Species Sensitivity from Target Interaction and Evolutionary Lineage**

SENTINEL models cross-species chemical susceptibility by tracing evolutionary divergence in xenobiotic target proteins. By integrating comparative sequence analysis, 3D structural modelling, and binding simulations, the workflow translates molecular initiating events (MIEs) into predictive Species Sensitivity Distributions (SSDs) across uncharacterised taxa.

The pipeline is engineered for SLURM-managed HPC clusters, distributing GPU-accelerated structural predictions and parallel CPU tasks across cluster nodes.

---

## Pipeline Architecture

```mermaid
flowchart TD
    %% Inputs
    subgraph INP["Input Data"]
        direction TB
        I1["Target Sequences / Accessions"]
        I2["Chemical Structures (.sdf / SMILES)"]
        I3["Empirical Ecotoxicity Endpoints (LC50/EC50)"]
    end

    %% Checkpoints
    subgraph EP["Execution Checkpoints (--entry)"]
        direction TB
        E_SEQ["--entry sequences (Default)"]
        E_STR["--entry structures"]
        E_POS["--entry poses"]
        E_ANA["--entry analysis"]
    end

    %% Stages
    subgraph ST1["Phase 1: Orthology & Lineage Inference (SLURM CPU)"]
        direction TB
        S1A["Target Homology & Domain Slicing<br/><b>HMMER v3.4</b>"]
        S1B["Multiple Sequence Alignment<br/><b>MAFFT v7.5</b>"]
        S1C["Active-Site Residue Calling & QC"]
        S1D{"Phylogeny Engine"}
        S1E["Lineage Tree<br/><b>FastTree</b>"]
        S1F["Bootstrapped ML Tree<br/><b>RAxML-NG</b> (--with-ml-tree)"]

        S1A --> S1B --> S1C --> S1D
        S1D -->|Default| S1E
        S1D -->|Optional| S1F
    end

    subgraph ST2["Phase 2: High-Throughput Structural Modelling (SLURM GPU)"]
        direction TB
        S2A["Ensemble Folding & Pocket Generation<br/><b>LocalColabFold (AlphaFold2)</b>"]
        S2B["Structural QC (pLDDT & PAE Filtering)"]
        S2C["Explicit-Solvent Minimisation<br/><b>GROMACS</b> (--with-minimise)"]

        S2A --> S2B
        S2B -.->|Optional| S2C
    end

    subgraph ST3["Phase 3: Molecular Docking & Binding Profiles (SLURM CPU)"]
        direction TB
        S3A["Ligand & Receptor Preparation<br/><b>RDKit / Meeko / Open Babel</b>"]
        S3B["Grid Definition & Auto-Box Geometry"]
        S3C["Conformational Docking Arrays<br/><b>AutoDock Vina</b>"]
        S3D["Pose Clustering & Affinity Extraction"]

        S3A --> S3C
        S3B --> S3C
        S3C --> S3D
    end

    subgraph ST4["Phase 4: Comparative Evolutionary Modelling (Master Node)"]
        direction TB
        S4A["Endpoint Harmonisation & Data Joining"]
        S4B["Phylogenetic Signal Estimation<br/><b>Pagel's &lambda; / Blomberg's K</b>"]
        S4C["Evolutionary Regressions & Sensitivity Models<br/><b>PGLS / Scikit-learn</b>"]
        S4D["Leave-One-Clade-Out Cross-Validation"]

        S4A --> S4B --> S4C --> S4D
    end

    %% Outputs
    subgraph OUT["Primary Outputs"]
        direction TB
        O1["master_table.csv<br/><i>(Alignments, calls & metadata)</i>"]
        O2["MODEL/univariate_results.tsv<br/><i>(Signal & clade error metrics)</i>"]
        O3["DASHBOARD/index.html<br/><i>(Interactive report & SSD curves)</i>"]
    end

    %% Wiring Inputs to Entry Points
    I1 --> E_SEQ
    E_SEQ --> S1A
    E_STR --> S2A
    E_POS --> S3A
    E_ANA --> S4A

    %% Inter-stage Links
    S1E --> S2A
    S1F --> S2A
    S2B --> S3A
    S2C --> S3A
    I2 --> S3A
    S3D --> S4A
    I3 --> S4A

    %% Final Outputs
    S4D --> O1
    S4D --> O2
    S4D --> O3

    %% Styling
    classDef io fill:#f8fafc,stroke:#64748b,stroke-width:1px;
    classDef slurm_cpu fill:#f0fdf4,stroke:#16a34a,stroke-width:1.5px;
    classDef slurm_gpu fill:#fef2f2,stroke:#dc2626,stroke-width:1.5px;
    classDef master fill:#eff6ff,stroke:#2563eb,stroke-width:1.5px;
    classDef checkpoint fill:#fefce8,stroke:#ca8a04,stroke-dasharray: 4 2;

    class I1,I2,I3,O1,O2,O3 io;
    class E_SEQ,E_STR,E_POS,E_ANA checkpoint;
    class S1A,S1B,S1C,S1D,S1E,S1F,S3A,S3B,S3C,S3D slurm_cpu;
    class S2A,S2B,S2C slurm_gpu;
    class S4A,S4B,S4C,S4D master;
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
