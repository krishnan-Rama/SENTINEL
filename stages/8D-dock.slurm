#!/bin/bash
#SBATCH --job-name=vina_dock
#SBATCH --partition=epyc
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/dock_%A_%a.out
#SBATCH --error=logs/dock_%A_%a.err

# Vina docking, one array task per receptor, all ligands inside.
#
#   sbatch --array=0-81 8D-dock.slurm RECEPTORS/pdbqt LIGANDS DOCKED
#
# Each receptor gets its OWN box, read from RECEPTORS/receptor_prep.tsv. The
# gorge axis differs between species, which is the entire point of the study; a
# shared box would impose one species' geometry on all 82 and manufacture
# exactly the kind of artefactual uniformity the analysis is meant to detect.
#
# --exhaustiveness 32 rather than the default 8. The gorge is a deep narrow
# channel and the default search under-samples it, which shows up as run to run
# scatter of several kJ/mol on the same input. --seed is fixed so the run is
# reproducible, and 3 independent seeds are run per pair so the docking noise
# can be quantified the same way the structural noise was. A score whose
# between-seed spread exceeds the between-species spread is not a measurement.

set -euo pipefail

module purge
module load micromamba/2.8.1
export PATH="${DOCK_ENV:-/mnt/ecotox/GROUP-smbpk/c23048124/envs/dock}/bin:$PATH"

REC_DIR=$(readlink -f "${1:?usage: sbatch --array=0-N 8D-dock.slurm <receptor_pdbqt_dir> <ligand_dir> <outdir>}")
LIG_DIR=$(readlink -f "${2:?}")
OUT=$(readlink -f "${3:?}")
PREP_TABLE=${PREP_TABLE:-$(dirname "$REC_DIR")/receptor_prep.tsv}
QC_TABLE=${QC_TABLE:-$REC_DIR/receptor_pdbqt_qc.tsv}
SEEDS=${SEEDS:-"1 2 3"}
EXH=${EXH:-32}
NMODES=${NMODES:-9}
CPUS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$OUT" logs
command -v vina >/dev/null || { echo "FATAL: vina not on PATH"; exit 1; }

# only receptors that passed the tautomer check
# NOTE: the TSVs are CRLF (Python csv.writer defaults to \r\n line endings).
# Python readers strip that transparently; awk does not, and a trailing \r on
# the last field silently corrupts it. Strip it explicitly everywhere.
mapfile -t READY < <(awk -F'\t' '{gsub(/\r/,"")} NR>1 && $2=="ok" {print $1}' \
                     "$QC_TABLE" | sort)
N=${#READY[@]}
TASK=${SLURM_ARRAY_TASK_ID:-0}
if [[ "$TASK" -ge "$N" ]]; then
    echo "task $TASK beyond $N ready receptors; nothing to do"; exit 0
fi
SP=${READY[$TASK]}
REC="$REC_DIR/$SP.pdbqt"
[[ -s "$REC" ]] || { echo "FATAL: $REC missing"; exit 1; }

# per-species box
read -r CX CY CZ SX SY SZ < <(awk -F'\t' -v s="$SP" '
  {gsub(/\r/,"")}
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
  $h["species"]==s { print $h["box_center_x"], $h["box_center_y"],
                           $h["box_center_z"], $h["box_size_x"],
                           $h["box_size_y"], $h["box_size_z"] }' "$PREP_TABLE")
[[ -n "${CX:-}" ]] || { echo "FATAL: no box for $SP in $PREP_TABLE"; exit 1; }

echo "=== $SP ==="
echo "  box centre ($CX, $CY, $CZ)  size ($SX, $SY, $SZ)"

SPOUT="$OUT/$SP"; mkdir -p "$SPOUT"
SUMMARY="$SPOUT/scores.tsv"
echo -e "species\tligand\tseed\tmode\taffinity_kcal_mol\trmsd_lb\trmsd_ub" > "$SUMMARY"

shopt -s nullglob
for LIG in "$LIG_DIR"/*.pdbqt; do
    LNAME=$(basename "$LIG" .pdbqt)
    for S in $SEEDS; do
        POSE="$SPOUT/${LNAME}_seed${S}.pdbqt"
        LOG="$SPOUT/${LNAME}_seed${S}.log"
        if [[ ! -s "$POSE" ]]; then
            vina --receptor "$REC" --ligand "$LIG" \
                 --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
                 --size_x "$SX" --size_y "$SY" --size_z "$SZ" \
                 --exhaustiveness "$EXH" --num_modes "$NMODES" \
                 --seed "$S" --cpu "$CPUS" \
                 --out "$POSE" > "$LOG" 2>&1 || {
                     echo "  FAILED $LNAME seed $S"; continue; }
        fi
        awk -v sp="$SP" -v lg="$LNAME" -v sd="$S" '
            /^[[:space:]]*[0-9]+[[:space:]]+-?[0-9]/ {
                print sp"\t"lg"\t"sd"\t"$1"\t"$2"\t"$3"\t"$4 }' "$LOG" \
            >> "$SUMMARY"
        BEST=$(awk '/^[[:space:]]*1[[:space:]]/ {print $2; exit}' "$LOG")
        echo "  $LNAME seed $S: best ${BEST:-NA} kcal/mol"
    done
done

echo
echo "scores: $SUMMARY"
echo
echo "Do not read the affinities as potencies. Vina scores the non-covalent"
echo "Michaelis complex; organophosphate potency is k2/Kd and nothing here"
echo "touches k2. The oxon-versus-parent and oxon-versus-TCP contrasts are the"
echo "informative comparisons, not the absolute numbers."
