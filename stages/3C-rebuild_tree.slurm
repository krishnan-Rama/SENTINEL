#!/bin/bash
#SBATCH --job-name=ache_phylo
#SBATCH --partition=epyc
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/phylo_%j.out
#SBATCH --error=logs/phylo_%j.err

# Rebuild the AChE tree for publication rather than for clade assignment.
#
# What was wrong with the first tree
# ----------------------------------
# It was built from hmmalign match states. Match-state extraction DELETES every
# insert column, and the profile came from only 14 anchor sequences, six of them
# vertebrate. Any region not represented in that small, vertebrate-leaning model
# became an insert and was discarded for all 276 sequences. The regions most
# likely to resolve relationships within Pancrustacea are precisely the ones
# that diverge from a vertebrate-dominated seed, so they were the first to go.
#
# This realigns the sequences themselves with L-INS-i, trims by gap fraction,
# and runs a proper ML search with bootstraps.
#
# What it will NOT fix
# --------------------
# A gene tree of a duplicated family groups by paralogue before taxonomy. Insect
# ace1 sequences will sit with other ace1 sequences rather than with the ace2 of
# the same species, and Hexapoda is nested inside Crustacea, so Malacostraca and
# Branchiopoda are not expected to be reciprocally monophyletic with Insecta.
# Intermixing of those three classes is the correct result, not a failure.
# If you need a tree for PGLS, use a dated SPECIES tree, not this.

set -euo pipefail

module purge
module load MAFFT/7.526-GCC-13.2.0-with-extensions
module load RAxML-NG/1.2.2-GCC-13.3.0
module load Python/3.11.5-GCCcore-13.2.0

ROOT=${ROOT:-$PWD}
CLASSIFY=${CLASSIFY:-$ROOT/CLASSIFY}
ANCHORS=${ANCHORS:-$ROOT/S1_sequences/anchors/reference_backbone.faa}
OUT=${OUT:-$ROOT/PHYLO}
THREADS=${SLURM_CPUS_PER_TASK:-16}
GAP_MAX=${GAP_MAX:-0.6}      # drop columns gapped in more than this fraction
BOOT=${BOOT:-200}
MODEL=${MODEL:-LG+G4+F}

mkdir -p "$OUT" "$ROOT/logs"

# ---- 1. collect the sequences, unaligned ----------------------------------
python3 - "$CLASSIFY" "$ANCHORS" "$OUT/seqs.faa" <<'PY'
import csv, os, sys
classify, anchors, out = sys.argv[1:4]

def read_fasta(p):
    recs, h, s = [], None, []
    for l in open(p):
        if l.startswith(">"):
            if h: recs.append((h, "".join(s)))
            h, s = l[1:].rstrip("\n"), []
        elif l.strip(): s.append(l.strip())
    if h: recs.append((h, "".join(s)))
    return recs

keep = set()
with open(os.path.join(classify, "active_site_calls.tsv"), newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["call"].startswith("AChE") or r["call"] == "BChE":
            keep.add(r["sequence_id"])

seqs = {h.split()[0]: s for h, s in read_fasta(
    os.path.join(classify, "candidates.faa"))}
seen, n = set(), 0
with open(out, "w") as fh:
    for h, s in read_fasta(anchors):
        tag = "ANCHOR_" + h.split("__")[0]
        if tag in seen: continue
        seen.add(tag)
        fh.write(f">{tag}\n{s}\n"); n += 1
    for k in sorted(keep):
        if k in seqs and k not in seen:
            seen.add(k)
            fh.write(f">{k}\n{seqs[k]}\n"); n += 1
print(f"  {n} unaligned sequences")
PY

# ---- 2. L-INS-i --------------------------------------------------------------
# --auto on this many sequences silently drops to FFT-NS-2, which is a
# progressive method with no iterative refinement. L-INS-i is the accurate
# setting and ~250 sequences of ~600 aa is well within its practical range.
if [[ ! -s "$OUT/aln.faa" ]]; then
    echo "aligning with L-INS-i (this is the slow step)"
    mafft --localpair --maxiterate 1000 --anysymbol \
          --thread "$THREADS" "$OUT/seqs.faa" > "$OUT/aln.faa"
fi

# ---- 3. trim gappy columns ---------------------------------------------------
python3 - "$OUT/aln.faa" "$OUT/aln.trim.faa" "$GAP_MAX" <<'PY'
import sys
inp, out, gmax = sys.argv[1], sys.argv[2], float(sys.argv[3])
recs, h, s = [], None, []
for l in open(inp):
    if l.startswith(">"):
        if h: recs.append((h, "".join(s)))
        h, s = l[1:].strip(), []
    elif l.strip(): s.append(l.strip())
if h: recs.append((h, "".join(s)))
n, L = len(recs), len(recs[0][1])
keep = [j for j in range(L)
        if sum(1 for _, q in recs if q[j] in "-.") / n <= gmax]
with open(out, "w") as fh:
    for hh, q in recs:
        t = "".join(q[j] for j in keep).upper().replace(".", "-")
        fh.write(f">{hh}\n")
        for i in range(0, len(t), 60): fh.write(t[i:i+60] + "\n")
print(f"  {n} sequences, {L} columns -> {len(keep)} kept "
      f"(gap fraction <= {gmax})")
PY

# ---- 4. ML tree with bootstraps ---------------------------------------------
cd "$OUT"
raxml-ng --check --msa aln.trim.faa --model "$MODEL" --prefix chk
raxml-ng --all --msa aln.trim.faa --model "$MODEL" \
         --prefix ache --threads "$THREADS" --seed 12345 \
         --bs-trees "$BOOT" --tree pars{10},rand{10}

echo
echo "best tree with supports: $OUT/ache.raxml.support"
echo "  model $MODEL, $BOOT bootstraps"
echo
echo "Root it on the CCE and non-catalytic anchors before interpreting. Nodes"
echo "below ~70 bootstrap should not be read as resolved relationships, and"
echo "within Pancrustacea a lot of them will be: that is a statement about how"
echo "much information one ~600 aa gene carries, not a fixable parameter."

