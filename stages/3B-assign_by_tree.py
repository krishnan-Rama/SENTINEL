#!/usr/bin/env python3
"""
Assign AChE candidates by tree topology, and pick one representative per species.

Why this exists
---------------
The residue classifier tests the acyl pocket in three states and I named those
states taxonomically: Phe/Phe as "vertebrate", non-Phe/Phe as "invertebrate".
That was wrong. The residues are a real observation; the labels assert a lineage
claim two positions cannot carry. The visible consequence in your output is six
teleosts carrying an `AChE_invertebrate` call, which is not a thing.

Clade membership is the evidence for what a sequence is. This script builds one
tree from every AChE call plus the verified anchors, and assigns each candidate
to its nearest anchor by patristic distance, reporting the margin to the runner
up so borderline cases are visible rather than silently decided.

It also selects one representative per species, preferring the candidate that
sits closest to the anchor appropriate for that species' lineage. Where the tree
and the residue call disagree, both are reported and the row is flagged. Do not
resolve those by deleting a column.

Outputs
-------
  tree_assignments.tsv   per sequence: nearest anchor, distance, margin,
                         residue call, and whether the two agree
  representatives.tsv    one row per species with the chosen sequence
  representatives.faa    the chosen sequences
  ache.tree              the Newick tree, for you to look at

  module load MAFFT/7.526-GCC-13.2.0-with-extensions \\
              FastTree/2.2-GCCcore-14.2.0 Python/3.11.5-GCCcore-13.2.0
  python3 3B-assign_by_tree.py --classify CLASSIFY \\
      --anchors S1_sequences/anchors/reference_backbone.faa \\
      --outdir TREE --threads 16
"""

import argparse, csv, os, re, shutil, subprocess, sys
from collections import defaultdict

# Which anchor a species of a given class should sit nearest, if the sequence
# really is that species' synaptic acetylcholinesterase. A mismatch is not
# automatically an error, but it is always worth looking at.
# Acceptable AChE clades per class, as SETS. The previous version had entries
# only for vertebrates and Insecta, so every mollusc, crustacean, annelid and
# branchiopod fell through to a fallback and was counted as failing an
# expectation that had never been defined for it.
#
# Note the real limitation this exposes: the anchor set contains exactly one
# non-insect invertebrate AChE (C. elegans ace-1). Molluscan and crustacean
# AChEs have no close reference, so their placement is weaker than the
# vertebrate or insect assignments and should be reported as such rather than
# presented at equal confidence.
VERT = {"AChE_vert"}
INV = {"AChE_inv", "AChE_insect", "AChE_vert"}
EXPECTED = {
    "Actinopteri": VERT, "Amphibia": VERT, "Mammalia": VERT, "Aves": VERT,
    "Lepidosauria": VERT, "Chondrichthyes": VERT, "Chondrostei": VERT,
    "Hyperoartia": VERT, "Myxini": VERT,
    "Insecta": {"AChE_insect", "AChE_inv"},
    "Malacostraca": INV, "Branchiopoda": INV, "Copepoda": INV,
    "Hexanauplia": INV, "Arachnida": INV, "Chilopoda": INV, "Diplopoda": INV,
    "Gastropoda": INV, "Bivalvia": INV, "Cephalopoda": INV, "Polyplacophora": INV,
    "Polychaeta": INV, "Clitellata": INV, "Hirudinea": INV,
    "Anthozoa": INV, "Hydrozoa": INV, "Scyphozoa": INV,
    "Rhabditophora": INV, "Trematoda": INV, "Cestoda": INV,
    "Chromadorea": INV, "Enoplea": INV, "Echinoidea": INV, "Asteroidea": INV,
    "Holothuroidea": INV, "Ascidiacea": INV,
}
ACHE_GROUPS = {"AChE_vert", "AChE_inv", "AChE_insect"}
ANCHOR_GROUP = {
    "AChE_Torpedo_californica": "AChE_vert", "AChE_Homo_sapiens": "AChE_vert",
    "AChE_Mus_musculus": "AChE_vert", "AChE_Danio_rerio": "AChE_vert",
    "AChE_Caenorhabditis_elegans": "AChE_inv",
    "Ace1_Anopheles_gambiae": "AChE_insect",
    "Ace2_Drosophila_melanogaster": "AChE_insect",
    "BChE_Homo_sapiens": "BChE", "BChE_Mus_musculus": "BChE",
    "CES1_Homo_sapiens": "CCE", "CES2_Homo_sapiens": "CCE",
    "JHE_Heliothis_virescens": "CCE",
    "Neuroligin1_Homo_sapiens": "dead", "Thyroglobulin_Homo_sapiens": "dead",
}


def run(cmd, out=None):
    with (open(out, "w") if out else subprocess.DEVNULL) as fh:
        p = subprocess.run(cmd, stdout=fh if out else subprocess.DEVNULL,
                           stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.exit(f"FAILED: {' '.join(cmd)}\n{p.stderr[:1500]}")


# Excluded from the tree by default. Thyroglobulin and neuroligin are
# catalytically dead and so distant that they contribute long-branch noise
# without informing the AChE / BChE / CCE boundary, which is the only boundary
# this tree needs to resolve. The classifier already screens dead enzymes on
# triad integrity.
EXCLUDE_FROM_TREE = {"Thyroglobulin_Homo_sapiens"}


def match_states(a2m_seq):
    """A2M: uppercase and '-' are match states; lowercase and '.' are inserts."""
    return "".join(c for c in a2m_seq if c.isupper() or c == "-")


def read_fasta(path):
    recs, h, s = [], None, []
    for l in open(path):
        if l.startswith(">"):
            if h: recs.append((h, "".join(s)))
            h, s = l[1:].rstrip("\n"), []
        elif l.strip():
            s.append(l.strip())
    if h: recs.append((h, "".join(s)))
    return recs


# --------------------------------------------------------------- newick -----

class Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self):
        self.name = ""
        self.length = 0.0
        self.children = []
        self.parent = None


def parse_newick(text):
    """Minimal parser. Handles FastTree output: unrooted, branch lengths,
    optional numeric support values on internal nodes."""
    text = text.strip().rstrip(";")
    i = 0

    def parse_node():
        nonlocal i
        n = Node()
        if text[i] == "(":
            i += 1
            while True:
                c = parse_node()
                c.parent = n
                n.children.append(c)
                if text[i] == ",":
                    i += 1
                    continue
                if text[i] == ")":
                    i += 1
                    break
                sys.exit(f"newick parse error at {i}: {text[i-20:i+20]!r}")
        j = i
        while i < len(text) and text[i] not in "(),:;":
            i += 1
        n.name = text[j:i].strip().strip("'\"")
        if i < len(text) and text[i] == ":":
            i += 1
            j = i
            while i < len(text) and (text[i].isdigit() or text[i] in ".eE+-"):
                i += 1
            try:
                n.length = float(text[j:i])
            except ValueError:
                n.length = 0.0
        return n

    return parse_node()


def leaves_and_depths(root):
    """leaf name -> (node, {ancestor node: cumulative distance to it})"""
    out = {}
    stack = [(root, [])]
    while stack:
        n, path = stack.pop()
        if not n.children:
            d, acc = {}, 0.0
            for anc in reversed(path):
                acc += anc[1]
                d[anc[0]] = acc
            d_self = {}
            acc = n.length
            cur = n.parent
            while cur is not None:
                d_self[cur] = acc
                acc += cur.length
                cur = cur.parent
            out[n.name] = d_self
        else:
            for c in n.children:
                stack.append((c, path + [(n, n.length)]))
    return out


def patristic(d1, d2):
    """Distance between two leaves given their ancestor-distance maps."""
    common = set(d1) & set(d2)
    if not common:
        return float("inf")
    return min(d1[a] + d2[a] for a in common)


def subtree_anchors(root, anchor_leaf_names):
    """node -> set of anchor leaf names beneath it, by post-order traversal."""
    order, stack = [], [root]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(n.children)
    out = {}
    for n in reversed(order):
        if not n.children:
            out[n] = {n.name} & anchor_leaf_names
        else:
            s = set()
            for c in n.children:
                s |= out[c]
            out[n] = s
    return out


def nearest_clade_anchors(leaf, sub):
    """Walk up from the leaf to the first ancestor whose subtree contains an
    anchor, and return (anchors_in_that_clade, clade_size_in_nodes).

    This is the question that matters: which anchors share the smallest clade
    with this sequence. Patristic distance answers a different question and is
    misled by rate variation, so a fast-evolving sequence on a long branch
    drifts toward whichever anchor is closest in summed branch length even when
    it is nowhere near it topologically."""
    cur = leaf.parent
    depth = 0
    while cur is not None:
        if sub.get(cur):
            return sub[cur], depth
        cur = cur.parent
        depth += 1
    return set(), depth


# ----------------------------------------------------------------- main -----

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classify", required=True)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--inventory", default=None,
                    help="ODB/species_inventory.tsv, for the class of each "
                         "species. Without it the expectation check is skipped.")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--include-bche", action="store_true",
                    help="also place BChE calls, which makes the AChE clade "
                         "boundary visible in the tree")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    for t in ("mafft", "FastTree"):
        if not shutil.which(t) and not shutil.which(t.lower()):
            sys.exit(f"{t} not on PATH.\n  module load "
                     "MAFFT/7.526-GCC-13.2.0-with-extensions "
                     "FastTree/2.2-GCCcore-14.2.0")
    fasttree = "FastTree" if shutil.which("FastTree") else "fasttree"

    calls = {}
    with open(os.path.join(a.classify, "active_site_calls.tsv"), newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            calls[r["sequence_id"]] = r
    keep = {k: v for k, v in calls.items()
            if v["call"].startswith("AChE")
            or (a.include_bche and v["call"] == "BChE")}
    print(f"placing {len(keep)} sequences ({len(calls)} classified in total)")

    seqs = {h.split()[0]: s
            for h, s in read_fasta(os.path.join(a.classify, "candidates.faa"))}
    missing = [k for k in keep if k not in seqs]
    if missing:
        print(f"  WARNING: {len(missing)} ids not found in candidates.faa")

    # Reuse the profile alignment the classifier already produced. Realigning
    # these sequences de novo with MAFFT is actively harmful: the anchor set
    # contains thyroglobulin (2,768 aa) and neuroligin (863 aa), which against
    # ~600 aa cholinesterases inject roughly 2,200 columns that are gaps for
    # every other sequence. FastTree on a mostly-gap alignment produces noise,
    # which is how sequences with intact catalytic triads ended up nearest to
    # thyroglobulin. The A2M match states are a fixed 638-column frame in which
    # the extra thyroglobulin residues fall into insert columns and disappear.
    a2m_path = os.path.join(a.classify, "aligned.a2m")
    aln = os.path.join(a.outdir, "for_tree.aln")
    tre = os.path.join(a.outdir, "ache.tree")
    anchor_names = []

    if os.path.exists(a2m_path):
        print(f"using profile alignment {a2m_path}")
        wanted = set(keep)
        n_written = 0
        # aligned.a2m still holds every duplicate pub_gene_id, because hmmalign
        # aligned the full input and the classifier only collapsed duplicates in
        # its own dictionary afterwards. FastTree refuses an alignment with
        # non-unique names, so deduplicate on the way out, keeping the first.
        seen, n_dup = set(), 0
        with open(aln, "w") as fh:
            for h, sq in read_fasta(a2m_path):
                k = h.split()[0]
                if k in seen:
                    n_dup += 1
                    continue
                seen.add(k)
                ms = match_states(sq)
                if k.startswith("ANCHOR__"):
                    tag = k.split("__")[1] if "__" in k[8:] else k[8:]
                    if tag in EXCLUDE_FROM_TREE:
                        continue
                    nm = "ANCHOR_" + re.sub(r"[^A-Za-z0-9_]", "_", tag)
                    anchor_names.append((nm, tag))
                    fh.write(f">{nm}\n{ms}\n")
                    n_written += 1
                elif k in wanted:
                    fh.write(f">{k}\n{ms}\n")
                    n_written += 1
        print(f"  {n_written} sequences, "
              f"{len(anchor_names)} anchors, fixed-width columns"
              + (f", {n_dup} duplicate ids dropped" if n_dup else ""))
    else:
        print("WARNING: no aligned.a2m found; falling back to de novo MAFFT. "
              "Long outgroup anchors will distort the tree.")
        inp = os.path.join(a.outdir, "for_tree.faa")
        with open(inp, "w") as fh:
            for h, sq in read_fasta(a.anchors):
                tag = h.split("__")[0]
                if tag in EXCLUDE_FROM_TREE:
                    continue
                nm = "ANCHOR_" + re.sub(r"[^A-Za-z0-9_]", "_", tag)
                anchor_names.append((nm, tag))
                fh.write(f">{nm}\n{sq}\n")
            for k in keep:
                if k in seqs:
                    fh.write(f">{k}\n{seqs[k]}\n")
        run(["mafft", "--auto", "--anysymbol", "--thread", str(a.threads), inp], aln)

    if not os.path.exists(tre):
        print("building tree...")
        run([fasttree, "-lg", aln], tre)
    root = parse_newick(open(tre).read())
    depths = leaves_and_depths(root)
    leaf_node = {}
    stack = [root]
    while stack:
        n = stack.pop()
        if not n.children:
            leaf_node[n.name] = n
        stack.extend(n.children)
    print(f"tree: {len(depths)} leaves")

    anchors = {nm: tag for nm, tag in anchor_names if nm in depths}
    absent = [tag for nm, tag in anchor_names if nm not in depths]
    if absent:
        print(f"  WARNING: anchors missing from tree: {', '.join(absent)}")
    if not anchors:
        sys.exit("no anchors in the tree; cannot assign")

    cls_of = {}
    if a.inventory and os.path.exists(a.inventory):
        with open(a.inventory, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                t = (r.get("binomial") or r.get("organism_name", "")).replace(" ", "_")
                cls_of.setdefault(t, r.get("class", ""))

    anchor_leafset = set(anchors)
    sub = subtree_anchors(root, anchor_leafset)

    rows, by_sp = [], defaultdict(list)
    for k in keep:
        if k not in depths:
            continue
        # topology first: which anchors share the smallest clade with this leaf
        clade, depth = nearest_clade_anchors(leaf_node[k], sub)
        clade_groups = sorted({ANCHOR_GROUP.get(anchors[nm], "?")
                               for nm in clade})
        # branch length only breaks ties WITHIN that clade
        if clade:
            ds = sorted(((patristic(depths[k], depths[nm]), anchors[nm])
                         for nm in clade), key=lambda x: x[0])
        else:
            ds = sorted(((patristic(depths[k], depths[nm]), tag)
                         for nm, tag in anchors.items()), key=lambda x: x[0])
        best_d, best_a = ds[0]
        second_d, second_a = ds[1] if len(ds) > 1 else (float("inf"), "")
        grp = (clade_groups[0] if len(clade_groups) == 1
               else ANCHOR_GROUP.get(best_a, "?"))
        sp = k.split("__")[0]
        cls = cls_of.get(sp, "")
        exp = EXPECTED.get(cls, set())
        rc = calls[k]["call"]
        # do the tree and the residue call tell the same story?
        agree = ((grp == "AChE_vert" and rc == "AChE_vertebrate")
                 or (grp in ("AChE_inv", "AChE_insect")
                     and rc == "AChE_invertebrate"))
        row = {"sequence_id": k, "species": sp, "class": cls,
               "residue_call": rc, "residue_conf": calls[k]["confidence"],
               "clade_anchors": ";".join(sorted(anchors[nm] for nm in clade)),
               "clade_groups": ";".join(clade_groups),
               "clade_unambiguous": len(clade_groups) == 1,
               "nodes_to_clade": depth,
               "nearest_anchor": best_a, "anchor_group": grp,
               "dist": round(best_d, 4),
               "second_anchor": second_a, "second_dist": round(second_d, 4),
               "margin": round(second_d - best_d, 4),
               "expected_groups": ";".join(sorted(exp)) if exp else "undefined",
               "in_ache_clade": grp in ACHE_GROUPS,
               "meets_expectation": (grp in exp) if exp else None,
               "tree_vs_residue": "agree" if agree else "DISAGREE",
               "paralogue": calls[k].get("paralogue", "")}
        rows.append(row)
        by_sp[sp].append(row)

    with open(os.path.join(a.outdir, "tree_assignments.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    from collections import Counter
    n_amb = sum(1 for r in rows if not r["clade_unambiguous"])
    print(f"\nclades containing more than one anchor group: {n_amb}/{len(rows)}")
    print("\n== clade group ==")
    for k, v in Counter(r["anchor_group"] for r in rows).most_common():
        print(f"  {k}: {v}")
    print("\n== tree versus residue call ==")
    for k, v in Counter(r["tree_vs_residue"] for r in rows).most_common():
        print(f"  {k}: {v}")
    print("\n== what the tree changed ==")
    demoted = [r for r in rows
               if r["residue_call"].startswith("AChE")
               and not r["in_ache_clade"]]
    promoted = [r for r in rows
                if r["residue_call"] == "BChE" and r["in_ache_clade"]]
    print(f"  residue-called AChE reassigned out of the AChE clades: "
          f"{len(demoted)}")
    for k, v in Counter(r["anchor_group"] for r in demoted).most_common():
        print(f"    -> {k}: {v}")
    print(f"  residue-called BChE that sit inside an AChE clade: {len(promoted)}")
    print("  These are corrections, not conflicts. A fish has BChE and CES "
          "genes; the acyl-pocket rule over-called them as AChE and the tree "
          "puts them back.")

    # per-species view: the only one that decides n
    print("\n== per species ==")
    n_ok = n_wrong = n_none = 0
    wrong = []
    for sp, hits in sorted(by_sp.items()):
        cls = cls_of.get(sp, "")
        exp = EXPECTED.get(cls, set())
        ache = [h for h in hits if h["in_ache_clade"]]
        if not ache:
            n_none += 1
            wrong.append((sp, cls, "no sequence in any AChE clade"))
        elif exp and not any(h["anchor_group"] in exp for h in ache):
            n_wrong += 1
            wrong.append((sp, cls, "AChE clade but wrong one: "
                          + ";".join(sorted({h["anchor_group"] for h in ache}))))
        else:
            n_ok += 1
    print(f"  species with an AChE in the expected clade: {n_ok}")
    print(f"  species with an AChE in an unexpected clade: {n_wrong}")
    print(f"  species with no AChE-clade sequence at all:  {n_none}")
    for sp, cls, why in wrong[:20]:
        print(f"    {sp:30s} {cls:14s} {why}")
    if len(wrong) > 20:
        print(f"    ... and {len(wrong)-20} more")

    # ---- representatives ---------------------------------------------------
    reps = []
    for sp, hits in sorted(by_sp.items()):
        cls = cls_of.get(sp, "")
        exp = EXPECTED.get(cls, set())
        pool = [h for h in hits if h["anchor_group"] == exp] if exp else []
        basis = "nearest to lineage-appropriate anchor"
        if not pool:
            pool = [h for h in hits
                    if h["anchor_group"] in ("AChE_vert", "AChE_inv", "AChE_insect")]
            basis = "nearest to any AChE anchor; lineage expectation not met"
        if not pool:
            pool, basis = hits, "no AChE clade placement; REVIEW"
        best = sorted(pool, key=lambda h: h["dist"])[0]
        reps.append({**best, "n_candidates": len(hits), "selection_basis": basis})

    with open(os.path.join(a.outdir, "representatives.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(reps[0]))
        w.writeheader(); w.writerows(reps)
    with open(os.path.join(a.outdir, "representatives.faa"), "w") as fh:
        for r in reps:
            fh.write(f">{r['sequence_id']}\n")
            s = seqs[r["sequence_id"]]
            for i in range(0, len(s), 60):
                fh.write(s[i:i+60] + "\n")

    print(f"\nrepresentatives: {len(reps)} species -> {a.outdir}/representatives.faa")
    for k, v in Counter(r["selection_basis"] for r in reps).most_common():
        print(f"  {k}: {v}")
    print("\nDISAGREE rows are where two residues and the whole alignment tell "
          "different stories. The tree is the better evidence, but a large "
          "count means the residue rule is miscalibrated for that lineage and "
          "you should say so rather than quietly dropping the column.")


if __name__ == "__main__":
    main()
