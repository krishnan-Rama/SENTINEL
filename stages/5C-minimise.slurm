#!/bin/bash
#SBATCH --job-name=ache_em
#SBATCH --partition=epyc
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/em_%j.out
#SBATCH --error=logs/em_%j.err

# Solvated energy minimisation of the ColabFold models.
#
#   ROOT=$PWD sbatch 5C-minimise.slurm FOLD/models EM
#
# Read this before running it
# ---------------------------
# ColabFold's --amber flag already performs a restrained minimisation, and for
# RIGID-RECEPTOR docking that is usually sufficient. This script exists for two
# other reasons: it removes the position restraints so side chains relax
# without bias toward the predicted rotamer, and it puts the protein in
# explicit solvent with counter-ions, which matters if you later run MD.
#
# If your next step is Vina or smina on a rigid receptor, you can skip this and
# lose very little. If you are going on to MD, you need it.
#
# The thing that actually affects your docking is NOT minimisation
# ----------------------------------------------------------------
# It is the protonation state of the catalytic histidine. In acetylcholinesterase
# His440 (Torpedo numbering) must be neutral with the proton on Nd1, so that Ne2
# faces Ser200 Og and can accept its proton. That is HID in Amber naming, HISD in
# GROMACS. pdb2gmx assigns histidine tautomers from local hydrogen-bond geometry
# and gets the catalytic histidine wrong often enough that it must be checked,
# not assumed. Get it wrong and the nucleophile is not activated, and every
# docking score in the study inherits the error.
#
# This script reports which tautomer was assigned for every histidine near the
# catalytic serine. It does not silently override it, because forcing a state
# without looking is how the error becomes invisible.

set -euo pipefail

module purge
module load GROMACS/2025.3-foss-2025b

IN=$(readlink -f "${1:?usage: sbatch 5C-minimise.slurm <model_dir> <outdir>}")
OUT=$(readlink -f "${2:?usage: sbatch 5C-minimise.slurm <model_dir> <outdir>}")
THREADS=${SLURM_CPUS_PER_TASK:-16}
FF=${FF:-amber99sb-ildn}
WATER=${WATER:-tip3p}
BOXPAD=${BOXPAD:-1.2}          # nm from protein to box edge

mkdir -p "$OUT" logs

cat > "$OUT/em.mdp" <<'EOF'
integrator      = steep
emtol           = 500.0
emstep          = 0.01
nsteps          = 50000
nstlist         = 10
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF

REPORT="$OUT/em_report.tsv"
echo -e "model\tstatus\tatoms\tEpot_kJ_mol\tmax_force\tcatalytic_his_note" > "$REPORT"

shopt -s nullglob
# rank_001 only: the top-ranked model per sequence
for PDB in "$IN"/**/*_relaxed_rank_001_*.pdb "$IN"/*_relaxed_rank_001_*.pdb; do
    [[ -e "$PDB" ]] || continue
    NAME=$(basename "$PDB" .pdb | sed 's/_relaxed_rank_001.*//')
    W="$OUT/work/$NAME"
    FINAL="$OUT/minimised/${NAME}.pdb"
    [[ -s "$FINAL" ]] && { echo "cached $NAME"; continue; }
    mkdir -p "$W" "$OUT/minimised"
    echo "=== $NAME ==="

    if ! gmx pdb2gmx -f "$PDB" -o "$W/proc.gro" -p "$W/topol.top" \
            -i "$W/posre.itp" -ff "$FF" -water "$WATER" -ignh \
            > "$W/pdb2gmx.log" 2>&1; then
        echo -e "$NAME\tpdb2gmx_failed\t\t\t\t" >> "$REPORT"
        echo "  pdb2gmx failed, see $W/pdb2gmx.log"
        continue
    fi

    HIS=$(grep -c "HISE\|HISD\|HISH" "$W/topol.top" 2>/dev/null || echo 0)
    HISNOTE=$(grep -oE "HIS[DEH]" "$W/topol.top" 2>/dev/null | sort | uniq -c \
              | tr '\n' ' ' || echo "none")

    gmx editconf -f "$W/proc.gro" -o "$W/box.gro" -c -d "$BOXPAD" -bt cubic \
        > "$W/editconf.log" 2>&1
    gmx solvate -cp "$W/box.gro" -cs spc216.gro -o "$W/solv.gro" \
        -p "$W/topol.top" > "$W/solvate.log" 2>&1
    gmx grompp -f "$OUT/em.mdp" -c "$W/solv.gro" -p "$W/topol.top" \
        -o "$W/ions.tpr" -maxwarn 2 > "$W/grompp_ions.log" 2>&1
    echo SOL | gmx genion -s "$W/ions.tpr" -o "$W/ions.gro" -p "$W/topol.top" \
        -pname NA -nname CL -neutral > "$W/genion.log" 2>&1

    gmx grompp -f "$OUT/em.mdp" -c "$W/ions.gro" -p "$W/topol.top" \
        -o "$W/em.tpr" -maxwarn 2 > "$W/grompp_em.log" 2>&1
    gmx mdrun -v -deffnm "$W/em" -ntomp "$THREADS" > "$W/mdrun.log" 2>&1

    EPOT=$(awk '/Potential Energy/ {print $4}' "$W/mdrun.log" | tail -1)
    FMAX=$(awk '/Maximum force/ {print $4}' "$W/mdrun.log" | tail -1)
    NAT=$(head -2 "$W/em.gro" | tail -1 | tr -d ' ')

    # strip solvent, keep protein only, for docking
    echo Protein | gmx trjconv -s "$W/em.tpr" -f "$W/em.gro" -o "$FINAL" \
        -pbc mol -ur compact > "$W/trjconv.log" 2>&1

    echo -e "$NAME\tok\t$NAT\t${EPOT:-NA}\t${FMAX:-NA}\t${HISNOTE}" >> "$REPORT"
    echo "  Epot ${EPOT:-NA} kJ/mol, Fmax ${FMAX:-NA}, His states: ${HISNOTE}"
done

echo
echo "minimised protein-only structures: $OUT/minimised/"
echo "report: $REPORT"
echo
echo "CHECK THE catalytic_his_note COLUMN. Every model needs the catalytic"
echo "histidine as HISD (neutral, proton on Nd1). If pdb2gmx assigned HISE"
echo "there, the Ne2 lone pair is not available to Ser Og, the nucleophile is"
echo "not activated, and the resulting docking geometry is wrong in the one"
echo "place this whole study depends on. Identify the catalytic His by its"
echo "distance to the catalytic Ser Og, not by residue number, since numbering"
echo "differs across species."
