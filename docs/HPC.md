# Running on a cluster

Two stages are expensive and both are parallel across species. The rest runs on a laptop.

## Structure prediction, stage 5B

```
bin/sentinel --gpu-partition gpu --time-fold 24:00:00 --max-parallel 4 --only 5B run
```

which resolves to an `sbatch --array=0-N` submission sized from the number of chunks. Use
`--hpc local` to print the command instead of submitting it.

One array task per chunk of sequences rather than one long job. Five models with templates
and six recycles will not reliably finish a full species set inside one wall-clock window,
and a single job that dies at hour 40 loses everything not yet written. `colabfold_batch`
skips sequences whose outputs already exist, so a timed-out task resumes on requeue.

Requests one GPU, 8 CPUs and 64 GB, with a 24 hour limit. Twelve sequences per chunk fits
that comfortably on an A100.

The one failure worth knowing about in advance: with `--templates`, MMseqs2 fetches template
structures from the server and ColabFold then runs `hhsearch` locally. If `hhsearch` is not
on PATH, every sequence dies at the template step, but only after the MSA search has
finished, so the job burns roughly 40 minutes per chunk and exits zero having produced
nothing. The script checks for the binary up front and refuses to start without it. Set
`TEMPLATES=0` to fold without templates, at a real cost in accuracy for a fold this well
covered in the PDB.

## Docking, stage 8D

```
bin/sentinel --partition epyc --time-dock 12:00:00 --only 8D run
```

The array is sized from the number of receptors with status `ok` in the tautomer check, so
there is no array bound to type by hand.

One task per receptor, all ligands inside, 8 CPUs and 8 GB each. Each receptor uses its own box from `receptor_prep.tsv`; a shared box would impose one species' gorge
geometry on all of them and manufacture exactly the artefactual uniformity the analysis
exists to detect.

## Phylogeny, stage 3C (optional)

```
bin/sentinel --with-ml-tree --only 3C --partition epyc run
```

L-INS-i alignment then RAxML-NG with 200 bootstraps, 16 CPUs and 32 GB, 12 hours. The
alignment is the slow step. This produces the publication tree; stage 3B does not wait for it.

## Minimisation, stage 5C (optional)

```
bin/sentinel --with-minimise --only 5C --partition epyc run
```

Explicit-solvent minimisation in GROMACS. ColabFold's Amber relaxation is already sufficient
for rigid-receptor docking, so skip this unless you are going on to molecular dynamics. If
you do run it, read the `catalytic_his_note` column of the report: `pdb2gmx` assigns
histidine tautomers from local hydrogen-bond geometry and gets the catalytic histidine wrong
often enough that it has to be checked rather than assumed.

## Adapting to another site

Partition, account, QOS and wall clock are all driver flags and override the `#SBATCH`
directives in the scripts, so those do not need editing. What does need editing is the
module load lines and the environment paths, which are from the reference system and will
not work anywhere else. If your site uses PBS or SGE rather than SLURM, run with
`--hpc local`: the driver prints the fully resolved command for each heavy stage and you
translate the submission wrapper around it.

The parts worth preserving if you rewrite a SLURM script are the array structure, the
resource requests, and the resume behaviour.

## Line endings

The TSVs are written by Python's csv module, which uses CRLF. Python readers strip that
transparently; `awk` does not, and a trailing carriage return on the last field corrupts it
silently. The SLURM scripts strip it explicitly. If you write your own shell tooling around
these tables, do the same.
