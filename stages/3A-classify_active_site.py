#!/usr/bin/env python3
"""
Active-site positional classifier.

Takes the 9,558 OrthoDB candidates and separates AChE from BChE, carboxylesterase
and the catalytically dead members, using residues rather than annotation.

Why not just filter on the description column
---------------------------------------------
2,128 of your entries have a description that matches nothing, and OrthoDB
descriptions in non-model invertebrates are inherited from automated transfer.
Filtering on the word "acetylcholinesterase" would keep mislabelled CES and
discard correctly-folded AChE that nobody named. The residues are the evidence.

How the coordinate frame is established
---------------------------------------
All sequences are aligned to a profile built from your verified anchor set, so
every sequence shares one column frame. The Torpedo anchor is then used to map
Torpedo mature numbering (Ser200, Phe288, Phe290, Trp279 and so on) onto those
columns.

The UniProt entry P04058 includes a signal peptide, so its numbering is offset
from the mature numbering used in the structural literature. This script does
NOT assume the offset. It locates the catalytic serine by the nucleophile elbow
motif, sets offset = (its UniProt position) - 200, prints the value, and then
validates the whole mapping by checking that the anchors give the residues they
are supposed to give. If human AChE does not come back as Phe/Phe at 288/290 and
human BChE as Leu/Val, the mapping is wrong and the script stops. Do not edit
that check out.

Classification
--------------
  triad broken                      -> dead        (neuroligin, gliotactin)
  288=F and 290=F                   -> AChE_vertebrate
  288 not F and 290=F               -> AChE_invertebrate
  288 in {L,I,V,M} and 290 in {V,I,L,A,G} -> BChE
  otherwise                         -> CCE_or_unclear

Every diagnostic residue is written out per sequence. The call is a convenience
column; the residues are the data. Inspect anything marked low confidence.

  python3 3A-classify_active_site.py \\
      --candidates ODB/by_species \\
      --anchors S1_sequences/anchors/reference_backbone.faa \\
      --outdir CLASSIFY --threads 16
"""

# Deliberately stdlib-only, and deliberately WITHOUT tempfile. Loading the
# EasyBuild HMMER module pulls OpenSSL/3 onto LD_LIBRARY_PATH, and the system
# Python then links against a libcrypto that lacks the symbol version it was
# built for. tempfile imports random, which imports hashlib, which needs
# libcrypto, so an unused import was enough to kill the script on module load.
# None of the imports below touch hashlib or ssl.
import argparse, csv, glob, os, re, shutil, subprocess, sys
from collections import defaultdict

# Torpedo californica AChE, MATURE numbering, as used in the structural
# literature and in 1EA5.
TORPEDO_POS = {
    "Tyr70":  "Y", "Asp72":  "D", "Trp84":  "W", "Gly118": "G", "Gly119": "G",
    "Tyr121": "Y", "Glu199": "E", "Ser200": "S", "Ala201": "A", "Trp279": "W",
    "Phe288": "F", "Phe290": "F", "Glu327": "E", "Phe330": "F", "Tyr334": "Y",
    "His440": "H",
}
TRIAD = ("Ser200", "Glu327", "His440")
GORGE_AROMATIC = ("Trp84", "Tyr121", "Trp279", "Phe288", "Phe290", "Phe330",
                  "Tyr334", "Tyr70")
ELBOW = re.compile(r"G[EQ]S[AG]G")          # nucleophile elbow, Ser is index 2

# Insect ace1 / ace2 references. Two ace genes are present in all insects
# except the Cyclorrhapha suborder of Diptera. In mosquitoes, aphids and
# Lepidoptera, ace1 is the synaptic organophosphate target; in Drosophila and
# other higher Diptera, ace1 is lost and ace2 is the target. Both are called
# AChE_invertebrate by the acyl-pocket test, so they must be separated here or
# the Insecta arm is meaningless.
ACE1_TAG = "Ace1_Anopheles_gambiae"
ACE2_TAG = "Ace2_Drosophila_melanogaster"
ACE_MARGIN = 5.0     # percentage points below which the call is not trusted


def run(cmd, **kw):
    p = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if p.returncode != 0:
        sys.exit(f"FAILED: {' '.join(cmd)}\n{p.stderr[:2000]}")
    return p.stdout


def read_fasta(path):
    recs, h, s = [], None, []
    with open(path) as fh:
        for l in fh:
            if l.startswith(">"):
                if h: recs.append((h, "".join(s)))
                h, s = l[1:].rstrip("\n"), []
            elif l.strip():
                s.append(l.strip())
    if h: recs.append((h, "".join(s)))
    return recs


def match_states(a2m_seq):
    """A2M: uppercase and '-' are match states, lowercase and '.' are inserts.
    Stripping inserts leaves a fixed-length string shared by all sequences."""
    return "".join(c for c in a2m_seq if c.isupper() or c == "-")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True,
                    help="directory of per-species .faa, or a single fasta")
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--manifest", default=None,
                    help="extraction_manifest.tsv. Defaults to the one beside "
                         "--candidates. Guards against stale files from a "
                         "previous extraction.")
    ap.add_argument("--torpedo-tag", default="AChE_Torpedo_californica")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    # Check inputs here rather than letting MAFFT fail on a missing file and
    # emit twenty lines of shell errors that hide the actual cause.
    if not os.path.exists(a.anchors):
        sys.exit(f"anchors not found: {a.anchors}\n"
                 "  (resolved to {})\n".format(os.path.abspath(a.anchors))
                 + "  find it with:\n"
                   "    find . -name reference_backbone.faa")
    if os.path.getsize(a.anchors) < 1000:
        sys.exit(f"anchors file is {os.path.getsize(a.anchors)} bytes. "
                 "Rerun the anchor builder.")
    if not os.path.exists(a.candidates):
        sys.exit(f"candidates not found: {a.candidates}")

    for tool in ("mafft", "hmmbuild", "hmmalign"):
        if not shutil.which(tool):
            sys.exit(f"{tool} not on PATH.\n"
                     "  module load HMMER/3.4-gompi-2023b "
                     "MAFFT/7.526-GCC-13.2.0-with-extensions")

    # ---- gather candidates ------------------------------------------------
    merged = os.path.join(a.outdir, "candidates.faa")
    n = 0
    with open(merged, "w") as out:
        if os.path.isfile(a.candidates):
            srcs = [a.candidates]
        else:
            found = sorted(glob.glob(os.path.join(a.candidates, "*.faa")))
            # A bare glob picks up files left behind by earlier runs with a
            # different species selection. The manifest is the record of what
            # the current extraction actually wrote.
            man = a.manifest or os.path.join(os.path.dirname(
                a.candidates.rstrip("/")), "extraction_manifest.tsv")
            srcs = found
            if os.path.exists(man):
                listed = []
                with open(man, newline="") as fh:
                    for r in csv.DictReader(fh, delimiter="\t"):
                        fp = r.get("fasta", "")
                        if fp and os.path.exists(fp):
                            listed.append(fp)
                if listed:
                    stale = len(found) - len(listed)
                    if stale > 0:
                        print(f"WARNING: {len(found)} .faa on disk but the "
                              f"manifest lists {len(listed)}. Using the "
                              f"manifest and ignoring {stale} stale file(s) "
                              "from an earlier extraction.")
                        print(f"  manifest: {man}")
                    srcs = sorted(listed)
        if not srcs:
            sys.exit(f"no .faa found under {a.candidates}")
        for fp in srcs:
            for h, s in read_fasta(fp):
                out.write(f">{h}\n{s}\n")
                n += 1
    print(f"candidates: {n} sequences from {len(srcs)} file(s)")

    # ---- profile from anchors --------------------------------------------
    aln = os.path.join(a.outdir, "anchors.aln")
    hmm = os.path.join(a.outdir, "anchors.hmm")
    if not os.path.exists(hmm):
        with open(aln, "w") as fh:
            fh.write(run(["mafft", "--auto", "--thread", str(a.threads),
                          a.anchors]))
        run(["hmmbuild", "--cpu", str(a.threads), "-n", "CHE", hmm, aln])
    print("profile built from anchors")

    # ---- align anchors + candidates in one frame -------------------------
    combined = os.path.join(a.outdir, "combined.faa")
    with open(combined, "w") as out:
        for h, s in read_fasta(a.anchors):
            out.write(f">ANCHOR__{h.split()[0]}\n{s}\n")
        for h, s in read_fasta(merged):
            out.write(f">{h.split()[0]}\n{s}\n")

    a2m = os.path.join(a.outdir, "aligned.a2m")
    print("aligning (hmmalign)...")
    with open(a2m, "w") as fh:
        p = subprocess.run(["hmmalign", "--amino", "--outformat", "A2M",
                            hmm, combined], stdout=fh, stderr=subprocess.PIPE,
                           text=True)
    if p.returncode != 0:
        sys.exit(f"hmmalign failed:\n{p.stderr[:2000]}")

    a2m_recs = read_fasta(a2m)
    raw, seqs, dup = {}, {}, 0
    for h, sq in a2m_recs:
        k = h.split()[0]
        if k in seqs:
            dup += 1
            continue
        raw[k] = sq
        seqs[k] = match_states(sq)
    if dup:
        print(f"NOTE: {dup} duplicate sequence ids collapsed. OrthoDB carries "
              "multiple proteins per gene for some entries, and they share a "
              "pub_gene_id. First occurrence kept.")
    L = {len(v) for v in seqs.values()}
    if len(L) != 1:
        sys.exit(f"match-state lengths differ: {sorted(L)[:5]}. "
                 "A2M parsing assumption is wrong; stop and inspect aligned.a2m")
    print(f"aligned: {len(seqs)} sequences, {L.pop()} match states")

    # ---- locate Torpedo and calibrate numbering --------------------------
    tkey = next((k for k in seqs if k.startswith("ANCHOR__")
                 and a.torpedo_tag in k), None)
    if tkey is None:
        sys.exit(f"Torpedo anchor '{a.torpedo_tag}' not found in the alignment")
    traw = next(s for h, s in read_fasta(a.anchors) if a.torpedo_tag in h)

    m = ELBOW.search(traw)
    if not m:
        sys.exit("nucleophile elbow motif not found in the Torpedo anchor")
    ser_uniprot = m.start() + 3               # 1-based position of the Ser
    offset = ser_uniprot - 200
    print(f"Torpedo catalytic Ser at UniProt position {ser_uniprot} "
          f"({traw[m.start():m.start()+5]}) -> mature-numbering offset {offset}")
    if not (0 <= offset <= 60):
        sys.exit(f"offset {offset} is implausible. Check the Torpedo anchor.")

    # UniProt residue number -> match-state column, from the aligned Torpedo.
    #
    # This MUST walk the full A2M string, not the match-state string. In A2M,
    # lowercase characters are residues assigned to insert states and '.' is
    # padding; match_states() strips both. Counting residues over the stripped
    # string therefore skips every inserted residue and shifts the numbering by
    # one for each. Torpedo has one such residue upstream of the catalytic
    # serine, which is exactly the off-by-one that tripped the Glu199 check.
    tal = seqs[tkey]
    res2col = {}
    r = col = 0
    for ch in raw[tkey]:
        if ch == "-":                 # deletion within a match state
            col += 1
        elif ch == ".":               # padding within an insert column
            continue
        elif ch.isupper():            # residue occupying a match state
            r += 1
            res2col[r] = col
            col += 1
        else:                         # lowercase: residue in an insert column
            r += 1
            res2col[r] = None         # no match column exists for it

    colmap = {}
    for label, expect in TORPEDO_POS.items():
        up = int(re.search(r"\d+", label).group()) + offset
        if up not in res2col:
            sys.exit(f"{label} (UniProt {up}) is beyond the Torpedo sequence "
                     f"(length {len(traw)}). Offset calibration is wrong.")
        col = res2col[up]
        if col is None:
            sys.exit(f"{label} (UniProt {up}) fell in an insert column, so it "
                     "has no shared coordinate. Rebuild the profile with more "
                     "anchors so this region becomes a match state.")
        got = tal[col]
        if got != expect:
            sys.exit(f"MAPPING FAILED at {label}: expected {expect}, "
                     f"Torpedo gives {got}. Numbering is wrong. Stop.")
        colmap[label] = col
    print("Torpedo mapping verified at all 16 positions")

    # ---- validate on the other anchors -----------------------------------
    def residues(key):
        s = seqs[key]
        return {lab: s[c] for lab, c in colmap.items()}

    checks = [("AChE_Homo_sapiens", {"Phe288": "F", "Phe290": "F", "Trp279": "W"}),
              ("BChE_Homo_sapiens", {"Phe288": "L", "Phe290": "V"})]
    print("\nanchor validation:")
    fails = []
    for tag, want in checks:
        k = next((x for x in seqs if x.startswith("ANCHOR__") and tag in x), None)
        if k is None:
            print(f"  {tag}: NOT PRESENT (skipped)")
            continue
        got = residues(k)
        bad = {p: (got[p], v) for p, v in want.items() if got[p] != v}
        print(f"  {tag}: " + ", ".join(f"{p}={got[p]}" for p in want)
              + ("  OK" if not bad else f"  MISMATCH {bad}"))
        if bad:
            fails.append((tag, bad))
    if fails:
        sys.exit("\nAnchor validation failed. The positional mapping does not "
                 "reproduce known active-site residues, so no call it makes can "
                 "be trusted. Fix this before proceeding.")

    # ---- ace1 / ace2 references -------------------------------------------
    def anchor_key(tag):
        return next((x for x in seqs if x.startswith("ANCHOR__") and tag in x),
                    None)

    k1, k2 = anchor_key(ACE1_TAG), anchor_key(ACE2_TAG)
    ace_ok = k1 is not None and k2 is not None
    if not ace_ok:
        print("\nWARNING: ace1 and/or ace2 reference missing from the anchor "
              "set. Invertebrate AChE will not be split by paralogue, and every "
              "insect in your dataset is then a coin flip between the true "
              "organophosphate target and the wrong gene.")
    else:
        s1, s2 = seqs[k1], seqs[k2]

        def pid(q, ref):
            """Identity over columns where both are non-gap, in the shared
            match-state frame. Crude, but transparent and uses no extra tools."""
            n = m = 0
            for a, b in zip(q, ref):
                if a != "-" and b != "-":
                    n += 1
                    m += (a == b)
            return 100.0 * m / n if n else 0.0

        def paralogue(q):
            i1, i2 = pid(q, s1), pid(q, s2)
            d = abs(i1 - i2)
            if d < ACE_MARGIN:
                return "ambiguous", round(i1, 1), round(i2, 1), round(d, 1)
            return ("ace1" if i1 > i2 else "ace2"), round(i1, 1), round(i2, 1), round(d, 1)

    # ---- classify ---------------------------------------------------------
    BCHE_288 = set("LIVM")
    BCHE_290 = set("VILAG")

    def call(r):
        triad_ok = all(r[p] == TORPEDO_POS[p] for p in TRIAD)
        if not triad_ok:
            miss = [p for p in TRIAD if r[p] != TORPEDO_POS[p]]
            return "dead_or_broken", "high" if len(miss) > 1 else "medium"
        f288, f290 = r["Phe288"], r["Phe290"]
        if f288 == "F" and f290 == "F":
            return "AChE_vertebrate", "high" if r["Trp279"] == "W" else "medium"
        if f288 != "F" and f290 == "F":
            return "AChE_invertebrate", "high" if r["Trp84"] == "W" else "medium"
        if f288 in BCHE_288 and f290 in BCHE_290:
            return "BChE", "high" if r["Trp279"] != "W" else "medium"
        return "CCE_or_unclear", "low"

    out = os.path.join(a.outdir, "active_site_calls.tsv")
    labels = list(TORPEDO_POS)
    per_species = defaultdict(list)
    tally = defaultdict(int)
    para_tally = defaultdict(int)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sequence_id", "species", "call", "confidence",
                    "paralogue", "pid_ace1", "pid_ace2", "pid_margin",
                    "triad_intact", "n_gorge_aromatic", "n_gaps_at_sites"]
                   + labels)
        for key, s in seqs.items():
            if key.startswith("ANCHOR__"):
                continue
            r = {lab: s[c] for lab, c in colmap.items()}
            cl, conf = call(r)
            gaps = sum(1 for v in r.values() if v == "-")
            if gaps >= 4:
                cl, conf = "too_incomplete", "low"
            arom = sum(1 for p in GORGE_AROMATIC if r[p] in "FWY")
            sp = key.split("__")[0]
            par, p1, p2, pm = "", "", "", ""
            if cl == "AChE_invertebrate" and ace_ok:
                par, p1, p2, pm = paralogue(s)
                para_tally[par] += 1
            w.writerow([key, sp, cl, conf, par, p1, p2, pm,
                        all(r[p] == TORPEDO_POS[p] for p in TRIAD),
                        arom, gaps] + [r[l] for l in labels])
            tally[cl] += 1
            if cl.startswith("AChE"):
                per_species[sp].append((key, cl, conf, arom, gaps, par))

    print("\n== calls ==")
    for k, v in sorted(tally.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v}")
    if para_tally:
        print("\n== invertebrate AChE paralogue ==")
        for k, v in sorted(para_tally.items(), key=lambda x: -x[1]):
            print(f"  {k}: {v}")

    # ---- one representative AChE per species ------------------------------
    rep = os.path.join(a.outdir, "ache_per_species.tsv")
    raw = {h.split()[0]: s for h, s in read_fasta(merged)}
    with open(rep, "w", newline="") as fh, \
         open(os.path.join(a.outdir, "ache_representatives.faa"), "w") as ff:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "n_ache_candidates", "paralogues_present",
                    "selected", "call", "confidence", "n_gorge_aromatic",
                    "status"])
        n_manual = 0
        for sp, hits in sorted(per_species.items()):
            pars = sorted({h[5] for h in hits if h[5]})
            # Never auto-pick between ace1 and ace2. Which one is the synaptic
            # organophosphate target depends on the lineage, not on sequence
            # quality, and no scoring function here knows that.
            if len([p for p in pars if p in ("ace1", "ace2")]) > 1:
                status = "MANUAL: ace1 and ace2 both present"
                n_manual += 1
                for h in sorted(hits, key=lambda x: (x[4], -x[3])):
                    w.writerow([sp, len(hits), ";".join(pars), h[0], h[1],
                                h[2], h[3], status])
                    if h[0] in raw:
                        ff.write(f">{h[0]}\n")
                        q = raw[h[0]]
                        for i in range(0, len(q), 60):
                            ff.write(q[i:i + 60] + "\n")
                continue
            best = sorted(hits, key=lambda x: (x[4], -x[3],
                                               0 if x[2] == "high" else 1))[0]
            status = "auto" if len(hits) == 1 else "auto: multiple candidates"
            w.writerow([sp, len(hits), ";".join(pars), best[0], best[1],
                        best[2], best[3], status])
            if best[0] in raw:
                ff.write(f">{best[0]}\n")
                q = raw[best[0]]
                for i in range(0, len(q), 60):
                    ff.write(q[i:i + 60] + "\n")
    print(f"\nspecies with at least one AChE call: {len(per_species)}")
    print(f"  {out}\n  {rep}\n  {a.outdir}/ache_representatives.faa")
    print(f"species needing a manual paralogue decision: {n_manual}")
    print("\nRows marked MANUAL have both ace1 and ace2 and BOTH sequences were "
          "written out. Which is the synaptic organophosphate target depends on "
          "the lineage: ace1 in mosquitoes, aphids and Lepidoptera; ace2 in "
          "Drosophila and other higher Diptera, where ace1 is lost. No scoring "
          "function in this script knows that, so it does not guess.")
    print("\npid_margin below 5 means the ace1/ace2 identity scores are too "
          "close to separate by identity alone. Those need a tree, not a "
          "threshold.")


if __name__ == "__main__":
    main()
