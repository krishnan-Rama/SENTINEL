#!/usr/bin/env python3
"""
Prepare sequences for structure prediction.

Why trimming is not optional
----------------------------
The OrthoDB sequences are full precursors. Every one carries an N-terminal
signal peptide, and most carry a C-terminal region (the tetramerisation /
WAT-PRAD domain in vertebrates, the GPI-anchor attachment signal in invertebrate
AChE1) that is cleaved or membrane-associated and is not part of the catalytic
domain.

Folding them intact does three specific kinds of damage to this study:

  1. The signal peptide is a hydrophobic helix with no native context. AlphaFold
     places it somewhere, and in a fraction of models it lies across the gorge
     mouth. Any gorge volume or docking box computed from those models is
     measuring an artefact.
  2. The C-terminal region is disordered in the monomer and inflates the
     bounding box, which shifts a box centred on the structure rather than on
     the catalytic serine.
  3. Both regions have low pLDDT and drag down the whole-model score, so a
     model with an excellent catalytic domain looks poor and one with a lucky
     terminus looks fine. Whole-model pLDDT stops being comparable across
     species, which is exactly the comparison you need.

How the boundaries are set
--------------------------
Not by a sequence-length heuristic and not by SignalP. The profile alignment you
already have defines a shared coordinate frame, so the cut points are the
columns corresponding to Torpedo mature residues --start and --end. Every
sequence is cut at the same structural position, which is the only way the
resulting models are comparable.

Defaults span the mature catalytic domain as crystallised in 1EA5. Residues in
insert columns inside that span are kept; only material outside it is removed.

  module purge
  module load Python/3.11.5-GCCcore-13.2.0
  python3 5A-prepare_for_folding.py --aln CLASSIFY/aligned.a2m \\
      --fasta FINAL_SET/endpoint_final.faa FINAL_SET/prediction_final.faa \\
      --outdir FOLD --chunk 12
"""

import argparse, csv, math, os, re, sys
from collections import OrderedDict

ELBOW = re.compile(r"G[EQ]S[AG]G")


def read_fasta(p):
    recs, h, s = [], None, []
    for l in open(p):
        if l.startswith(">"):
            if h: recs.append((h, "".join(s)))
            h, s = l[1:].rstrip("\n"), []
        elif l.strip(): s.append(l.strip())
    if h: recs.append((h, "".join(s)))
    return recs


def match_states(q):
    return "".join(c for c in q if c.isupper() or c == "-")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True, help="CLASSIFY/aligned.a2m")
    ap.add_argument("--fasta", nargs="+", required=True,
                    help="FASTA or TSV listing the sequences to trim. TSVs "
                         "with a sequence_id column are preferred and exact; "
                         "FASTA headers are matched on the species tag, since "
                         "build_final_set writes them as Species|tier|clade "
                         "rather than as OrthoDB sequence ids.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--start", type=int, default=1,
                    help="first Torpedo MATURE residue to keep (default 1)")
    ap.add_argument("--end", type=int, default=537,
                    help="last Torpedo MATURE residue to keep (default 537, "
                         "the C-terminus of the 1EA5 catalytic domain)")
    ap.add_argument("--min-len", type=int, default=400)
    ap.add_argument("--max-len", type=int, default=700,
                    help="flag sequences still longer than this AFTER trimming. "
                         "The catalytic domain is ~537 residues; much more "
                         "means a large insertion inside the domain span, "
                         "usually a fused or mispredicted model.")
    ap.add_argument("--allow-tag-match", action="store_true",
                    help="permit matching on species tag when no exact "
                         "sequence ids are supplied. Off by default because it "
                         "selects every family member for each species.")
    ap.add_argument("--chunk", type=int, default=0,
                    help="split output into chunks of N sequences for array "
                         "jobs; 0 writes a single file")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    aln = OrderedDict()
    for h, s in read_fasta(a.aln):
        aln.setdefault(h.split()[0], s)

    tor = next((k for k in aln if "Torpedo" in k), None)
    if tor is None:
        sys.exit("no Torpedo anchor in the alignment; cannot set cut points")

    # residue -> match-state column, walking the FULL a2m so insert residues
    # are counted (they are residues; they just have no shared column)
    res2col, r, col = {}, 0, 0
    for ch in aln[tor]:
        if ch == "-":
            col += 1
        elif ch == ".":
            continue
        elif ch.isupper():
            r += 1; res2col[r] = col; col += 1
        else:
            r += 1; res2col[r] = None

    raw_tor = "".join(c for c in aln[tor] if c not in "-.").upper()
    m = ELBOW.search(raw_tor)
    if not m:
        sys.exit("nucleophile elbow not found in the Torpedo anchor")
    offset = (m.start() + 3) - 200
    print(f"Torpedo catalytic Ser at position {m.start()+3}; "
          f"mature-numbering offset {offset}")

    def col_of(mature):
        up = mature + offset
        c = res2col.get(up)
        while c is None and up in res2col:      # nudge past an insert residue
            up += 1
            c = res2col.get(up)
        return c

    c_start, c_end = col_of(a.start), col_of(a.end)
    if c_start is None or c_end is None or c_end <= c_start:
        sys.exit(f"could not resolve cut columns for mature {a.start}-{a.end}")
    print(f"keeping alignment columns {c_start} to {c_end} "
          f"(Torpedo mature {a.start}-{a.end})")

    # Two ways of naming the same sequences. TSVs from build_final_set carry
    # the OrthoDB sequence_id exactly; the matching FASTA headers were rewritten
    # to Species|tier|clade and carry only the species tag. Accept both, and
    # match a candidate if either its full id or its species tag was requested.
    want_ids, want_tags = set(), set()
    for fp in a.fasta:
        if fp.lower().endswith((".tsv", ".csv")):
            sep = "\t" if fp.lower().endswith(".tsv") else ","
            with open(fp, newline="") as fh:
                for r in csv.DictReader(fh, delimiter=sep):
                    sid = (r.get("sequence_id") or "").strip()
                    if sid:
                        want_ids.add(sid)
                        want_tags.add(sid.split("__")[0])
        else:
            for h, _ in read_fasta(fp):
                tok = h.split()[0]
                if "__" in tok:
                    want_ids.add(tok)
                    want_tags.add(tok.split("__")[0])
                else:
                    want_tags.add(tok.split("|")[0])
    if not (want_ids or want_tags):
        sys.exit("no sequence identifiers found in " + ", ".join(a.fasta))
    print(f"requested: {len(want_ids)} exact ids, {len(want_tags)} species tags")
    if not want_ids and not a.allow_tag_match:
        sys.exit(
            "\nSTOP: no exact sequence ids were found, only species tags.\n"
            "Matching on species tag selects EVERY cholinesterase-family member\n"
            "for each species, not the AChE representative you chose. On this\n"
            "dataset that is roughly 650 sequences instead of 85, and most of\n"
            "them are BChE, carboxylesterases and OrthoDB gene-id collisions.\n\n"
            "Use the TSVs, which carry the sequence_id column:\n"
            "  --fasta ../FINAL_SET/endpoint_final.tsv "
            "../FINAL_SET/prediction_final.tsv\n\n"
            "Pass --allow-tag-match only if you really do want every family\n"
            "member folded.")

    out, report, seen_seq = [], [], {}
    for key, a2m in aln.items():
        if key.startswith("ANCHOR__"):
            continue
        # Exact ids take PRECEDENCE. want_tags is also populated from the ids
        # (for reporting), so an "either matches" test lets the tag fallback
        # fire even when exact ids were supplied, which selects every family
        # member for each species instead of the chosen representative.
        if want_ids:
            if key not in want_ids:
                continue
        elif key.split("__")[0] not in want_tags:
            continue
        kept, cur, n_pre, n_post = [], 0, 0, 0
        for ch in a2m:
            if ch == ".":
                continue
            if ch == "-":
                cur += 1
                continue
            if ch.isupper():
                inside = c_start <= cur <= c_end
                if inside: kept.append(ch)
                elif cur < c_start: n_pre += 1
                else: n_post += 1
                cur += 1
            else:                                # insert residue
                if c_start <= cur <= c_end: kept.append(ch.upper())
                elif cur < c_start: n_pre += 1
                else: n_post += 1
        seq = "".join(kept)
        full = len(seq) + n_pre + n_post
        if len(seq) < a.min_len:
            flag = "SHORT"
        elif len(seq) > a.max_len:
            # large insertion inside the catalytic domain span: a fused model
            # or a mispredicted gene. Folding it wastes GPU time and the extra
            # mass sits somewhere unpredictable relative to the gorge.
            flag = "LONG"
        else:
            flag = "ok"
        report.append({"sequence_id": key, "full_length": full,
                       "trimmed_length": len(seq), "removed_N": n_pre,
                       "removed_C": n_post, "status": flag})
        if flag != "ok":
            continue
        if seq in seen_seq:                      # identical after trimming
            report[-1]["status"] = f"duplicate_of:{seen_seq[seq]}"
            continue
        seen_seq[seq] = key
        out.append((key, seq))

    if not report:
        ex_aln = [k for k in aln if not k.startswith("ANCHOR__")][:3]
        ex_want = sorted(want_ids)[:3] or sorted(want_tags)[:3]
        sys.exit("no sequences matched.\n"
                 f"  alignment keys look like: {ex_aln}\n"
                 f"  requested look like:      {ex_want}\n"
                 "  Pass the TSVs instead of the FASTAs:\n"
                 "    --fasta FINAL_SET/endpoint_final.tsv "
                 "FINAL_SET/prediction_final.tsv")

    with open(os.path.join(a.outdir, "trim_report.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(report[0]))
        w.writeheader(); w.writerows(report)

    def write(path, recs):
        with open(path, "w") as fh:
            for k, s in recs:
                fh.write(f">{k}\n")
                for i in range(0, len(s), 60):
                    fh.write(s[i:i+60] + "\n")

    if a.chunk and a.chunk > 0:
        d = os.path.join(a.outdir, "chunks"); os.makedirs(d, exist_ok=True)
        n = math.ceil(len(out) / a.chunk)
        for i in range(n):
            write(os.path.join(d, f"chunk_{i:03d}.fasta"),
                  out[i*a.chunk:(i+1)*a.chunk])
        print(f"\n{len(out)} sequences -> {n} chunks in {d}/")
        print(f"  submit with: sbatch --array=0-{n-1} 5B-fold.slurm {d} FOLD/models")
    else:
        write(os.path.join(a.outdir, "to_fold.fasta"), out)
        print(f"\n{len(out)} sequences -> {a.outdir}/to_fold.fasta")

    lens = [len(s) for _, s in out]
    nshort = sum(1 for r in report if r["status"] == "SHORT")
    nlong = sum(1 for r in report if r["status"] == "LONG")
    ndup = sum(1 for r in report if str(r["status"]).startswith("duplicate"))
    print(f"  trimmed length: min {min(lens)}, median "
          f"{sorted(lens)[len(lens)//2]}, max {max(lens)}")
    print(f"  removed N-terminal: median "
          f"{sorted(r['removed_N'] for r in report)[len(report)//2]} residues")
    print(f"  removed C-terminal: median "
          f"{sorted(r['removed_C'] for r in report)[len(report)//2]} residues")
    print(f"  excluded: {nshort} too short (<{a.min_len}), "
          f"{nlong} too long (>{a.max_len}), "
          f"{ndup} identical after trimming")
    print(f"\ntrim_report.tsv records what was cut from each sequence. If the "
          f"median N-terminal cut is not roughly 20-30 residues, check the "
          f"offset before folding anything.")


if __name__ == "__main__":
    main()
