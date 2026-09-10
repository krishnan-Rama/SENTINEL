#!/bin/bash
#SBATCH --job-name=colabfold
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=logs/colabfold_%A_%a.out
#SBATCH --error=logs/colabfold_%A_%a.err

# ColabFold for the AChE comparative set.
#
#   sbatch --array=0-6 5B-fold.slurm FOLD/chunks FOLD/models
#
# Array over chunks rather than one long job: 85 models with templates and five
# seeds will not reliably finish inside one wall-clock window, and a single job
# that dies at hour 40 loses everything not yet written.
#
# Differences from a default run, and why each one matters here
# -------------------------------------------------------------
# --templates       Cholinesterase is one of the best-covered folds in the PDB
#                   (1EA5, 4EY7, 1QO9 and many more). Withholding templates for
#                   no reason costs accuracy in exactly the region you measure.
# --num-models 5    You are comparing gorge geometry ACROSS species. A single
#                   seed gives one sample of side-chain placement, and the
#                   between-species differences you are looking for are of the
#                   same size as between-seed scatter. Five models per species
#                   let you show that the differences exceed the noise, or that
#                   they do not, which is itself a result.
# --num-recycle 6   Cheap, and helps convergence on the long gorge loops.
# --amber           Restrained minimisation. Necessary before docking: raw
#                   AlphaFold output routinely has minor clashes and
#                   non-physical bond geometry that a docking scoring function
#                   will happily interpret as structure.
#
# NOT set: --msa-mode single_sequence. Deep MSAs are the point for a
# comparative set; single-sequence mode would make the models depend on
# per-species alignment depth in an uncontrolled way.

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: sbatch --array=0-N 5B-fold.slurm <chunk_dir_or_fasta> <outdir>"
    exit 1
fi

IN=$(readlink -f "$1")
OUT=$(readlink -f "$2")
mkdir -p "$OUT" logs

PROJECT_ROOT="/mnt/ecotox/GROUP-smbpk/c23048124/colabfold"
PIXI_ENV="${PROJECT_ROOT}/localcolabfold/.pixi/envs/default"
COLABFOLD_BIN="${PIXI_ENV}/bin/colabfold_batch"
export LD_LIBRARY_PATH="${PIXI_ENV}/lib:${LD_LIBRARY_PATH:-}"
export XDG_CACHE_HOME="${PROJECT_ROOT}/cache"

[[ -x "$COLABFOLD_BIN" ]] || { echo "FATAL: $COLABFOLD_BIN not executable"; exit 1; }

# Put the pixi env's bin on PATH so its helper binaries are visible.
export PATH="${PIXI_ENV}/bin:${PATH}"

# --templates needs hhsearch. MMseqs2 fetches the template CIFs from the
# server, then ColabFold runs hhsearch LOCALLY to align and rank them. Without
# the binary every sequence dies at get_msa_and_templates, but only AFTER the
# MSA search has completed, so the job burns ~40 min per chunk and exits 0 with
# no structures. Check it up front.
TEMPLATES=${TEMPLATES:-1}
TEMPLATE_FLAG="--templates"
if [[ "$TEMPLATES" == "1" ]]; then
    if ! command -v hhsearch >/dev/null 2>&1; then
        echo "FATAL: --templates requested but hhsearch is not on PATH."
        echo "  install it into the pixi env:"
        echo "    module load micromamba/2.8.1"
        echo "    micromamba install -y -c bioconda -c conda-forge \\"
        echo "      -p ${PIXI_ENV} hhsuite"
        echo "  or rerun this job with TEMPLATES=0 to fold without templates."
        exit 1
    fi
    echo "hhsearch: $(command -v hhsearch)"
else
    TEMPLATE_FLAG=""
    echo "NOTE: templates disabled. Cholinesterase is exceptionally well"
    echo "covered in the PDB, so this costs real accuracy in the gorge."
fi

# GPU relaxation needs a CUDA-capable OpenMM. If it is missing, --amber still
# works on CPU; slower, but it does not silently produce nothing.
RELAX_FLAG="--amber --use-gpu-relax"
if [[ "${GPU_RELAX:-1}" != "1" ]]; then
    RELAX_FLAG="--amber"
fi

# pick this array task's chunk
if [[ -d "$IN" ]]; then
    TASK=${SLURM_ARRAY_TASK_ID:-0}
    FASTA=$(printf '%s/chunk_%03d.fasta' "$IN" "$TASK")
    SUB="$OUT/chunk_$(printf '%03d' "$TASK")"
else
    FASTA="$IN"; SUB="$OUT"
fi
[[ -s "$FASTA" ]] || { echo "FATAL: no input at $FASTA"; exit 1; }
mkdir -p "$SUB"

echo "input:  $FASTA  ($(grep -c '^>' "$FASTA") sequences)"
echo "output: $SUB"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || \
    echo "WARNING: no GPU visible; this will be extremely slow"

# colabfold_batch skips sequences whose outputs already exist, so requeuing a
# timed-out task resumes rather than restarting.
"${COLABFOLD_BIN}" \
    ${TEMPLATE_FLAG} \
    --num-models 5 \
    --num-recycle 6 \
    --rank plddt \
    ${RELAX_FLAG} \
    --model-type auto \
    "$FASTA" "$SUB"

N_PDB=$(ls "$SUB"/*_relaxed_rank_001_*.pdb 2>/dev/null | wc -l)
N_SEQ=$(grep -c '^>' "$FASTA")
echo "done: $N_PDB relaxed rank-1 models from $N_SEQ sequences"
if [[ "$N_PDB" -lt "$N_SEQ" ]]; then
    echo "WARNING: fewer models than sequences. Check $SUB/log.txt for the"
    echo "reason; a chunk can exit 0 having produced nothing at all."
fi
echo
echo "Next: QC before anything else. A model with a well-folded scaffold and a"
echo "disordered acyl pocket is worse than useless here, because it will dock"
echo "and give a number. Check per-residue pLDDT at the 16 active-site"
echo "positions, not the whole-model mean."
