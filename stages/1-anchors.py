#!/usr/bin/env python3
"""
S1a: Rebuild the reference backbone by verified accession.

Replaces the gene-name query approach in the first version of S1, which
silently returned P02795 (metallothionein-2, 61 aa) for the CES1 slot because
UniProt's gene field matched loosely. That sequence would have been written
into the backbone and used to root the AChE / BChE / CCE placement in S3.

Two changes:
  1. Anchors are fetched by accession, not by gene-name query.
  2. Every anchor is checked against an expected length window and a required
     substring of the protein name. A mismatch is a hard failure, not a warning.

Accessions below were confirmed from your own S1 run except CES1_HUMAN, which
was wrong and is corrected to P23141.

Usage:
  python3 1-anchors.py --outdir S1_sequences --contact you@cardiff.ac.uk
"""

import argparse
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

UNIPROT = "https://rest.uniprot.org"
FIELDS = "accession,reviewed,id,protein_name,gene_names,organism_name,length,sequence"

# tag: (accession, min_len, max_len, required substring in protein name, role)
# role drives how S3 uses the sequence: 'target' clade, 'outgroup' clade, or
# 'dead' (catalytically inactive, defines the floor of the family).
ANCHORS = {
    "AChE_Torpedo_californica":     ("P04058", 560, 610, "cholinesterase", "AChE"),
    "AChE_Homo_sapiens":            ("P22303", 590, 640, "cholinesterase", "AChE"),
    "AChE_Mus_musculus":            ("P21836", 590, 640, "cholinesterase", "AChE"),
    "AChE_Danio_rerio":             ("Q9DDE3", 600, 660, "cholinesterase", "AChE"),
    "AChE_Caenorhabditis_elegans":  ("P38433", 590, 650, "cholinesterase", "AChE"),
    "Ace2_Drosophila_melanogaster": ("P07140", 620, 680, "cholinesterase", "AChE_ace2"),
    "Ace1_Anopheles_gambiae":       ("Q869C3", 690, 780, "cholinesterase", "AChE_ace1"),
    "BChE_Homo_sapiens":            ("P06276", 580, 630, "cholinesterase", "BChE"),
    "BChE_Mus_musculus":            ("Q03311", 580, 630, "cholinesterase", "BChE"),
    "CES1_Homo_sapiens":            ("P23141", 540, 590, "esterase",       "CCE"),
    "CES2_Homo_sapiens":            ("O00748", 530, 580, "esterase",       "CCE"),
    "JHE_Heliothis_virescens":      ("P12992", 540, 590, "esterase",       "CCE"),
    "Neuroligin1_Homo_sapiens":     ("Q8N2Q7", 800, 900, "neuroligin",     "dead"),
    "Thyroglobulin_Homo_sapiens":   ("P01266", 2600, 2900, "thyroglobulin", "dead"),
}

# Apis mellifera was returned as A0ACM8Q821 (TrEMBL, 628 aa) in your run. Bees
# carry both ace1 and ace2 and the TrEMBL name does not say which. It is
# excluded here deliberately. Add it back only after S4 assigns the paralogue.


def http_get(url, contact, tries=5, timeout=90):
    hdr = {"User-Agent": f"ache-era-pipeline/1.1 ({contact})"}
    delay = 2.0
    for _ in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=hdr),
                                        timeout=timeout) as r:
                return r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return None
            time.sleep(float(e.headers.get("Retry-After", delay)))
        except Exception:
            time.sleep(delay)
        delay *= 2
    return None


def fetch(acc, contact):
    url = (f"{UNIPROT}/uniprotkb/search?"
           + urllib.parse.urlencode({"query": f"accession:{acc}", "format": "tsv",
                                     "fields": FIELDS, "size": 1}))
    txt = http_get(url, contact)
    if not txt:
        return None
    lines = txt.rstrip("\n").split("\n")
    if len(lines) < 2:
        return None
    return dict(zip(lines[0].split("\t"), lines[1].split("\t")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--contact", required=True)
    a = ap.parse_args()

    d = os.path.join(a.outdir, "anchors")
    os.makedirs(d, exist_ok=True)
    fa = os.path.join(d, "reference_backbone.faa")
    tb = os.path.join(d, "reference_backbone.tsv")

    good, bad = [], []
    with open(tb, "w") as t:
        t.write("tag\trole\taccession\tuniprot_id\tprotein_name\torganism\tlength\treviewed\tverdict\n")
        for tag, (acc, lo, hi, need, role) in ANCHORS.items():
            r = fetch(acc, a.contact)
            if r is None:
                print(f"  FAIL {tag}: accession {acc} not retrieved")
                bad.append(tag)
                continue
            L = int(r.get("Length", "0") or 0)
            pname = r.get("Protein names", "")
            ok_len = lo <= L <= hi
            ok_name = need.lower() in pname.lower()
            verdict = "pass" if (ok_len and ok_name) else \
                      ("FAIL_length" if not ok_len else "FAIL_name")
            t.write(f"{tag}\t{role}\t{acc}\t{r.get('Entry Name','')}\t{pname}\t"
                    f"{r.get('Organism','')}\t{L}\t{r.get('Reviewed','')}\t{verdict}\n")
            if verdict == "pass":
                good.append((tag, role, r))
                print(f"  ok   {tag}: {acc} {L} aa  {pname[:55]}")
            else:
                bad.append(tag)
                print(f"  FAIL {tag}: {acc} {L} aa ({verdict}) {pname[:55]}")
            time.sleep(0.35)

    with open(fa, "w") as f:
        for tag, role, r in good:
            f.write(f">{tag}__{r['Entry']}__{role}\n")
            s = r.get("Sequence", "")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")

    print(f"\n{len(good)}/{len(ANCHORS)} anchors verified -> {fa}")
    print(f"verification table -> {tb}")
    if bad:
        print("FAILED: " + ", ".join(bad))
        print("Do not run S3 until these are resolved. An unverified backbone "
              "makes every downstream clade assignment unfalsifiable.")
        sys.exit(1)


if __name__ == "__main__":
    main()
