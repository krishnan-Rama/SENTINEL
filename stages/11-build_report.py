#!/usr/bin/env python3
"""
Build a single-file tiered dashboard from the pipeline outputs.

Design constraint that shapes everything else
---------------------------------------------
A dashboard that presents a species ranking without the evidence for whether
that ranking is trustworthy would automate the failure this pipeline exists to
detect. So the dashboard is GATED: the model-validation status is computed
first, and the prediction tab states plainly which predictors survived
correction and which did not. If none survived, it says so at the top, in the
banner, before any table.

What it predicts, and why that is not the structural descriptors
----------------------------------------------------------------
When phylogenetic signal in the endpoint is high, the best available predictor
of an unmeasured species is its relatives, not its protein geometry. Prediction
here is therefore phylogenetic imputation under Brownian motion with Pagel's
lambda: the conditional expectation of an unobserved tip given the observed
tips and the tree, with a credible interval from the conditional variance. That
is standard read-across, made quantitative. It is honest about what the data
supports and it is useful to a regulator in a way a non-validated structural
ranking is not.

Output is one HTML file with the data embedded. No server, no CDN, no network.
It can be emailed, archived alongside a dossier, or opened offline.

  module purge
  module load Python/3.11.5-GCCcore-13.2.0 SciPy-bundle/2025.07-gfbf-2025b
  python3 11-build_report.py --root . --outdir DASHBOARD
"""

import argparse, csv, html, json, math, os, re, sys, time
from collections import defaultdict, Counter

try:
    import numpy as np
except ImportError:
    sys.exit("numpy required. module load SciPy-bundle/2025.07-gfbf-2025b")


# Functional annotation of the diagnostic positions, Torpedo mature numbering,
# with the mammalian equivalent in brackets. These are what make the active-site
# string in the species panel interpretable rather than decorative.
RESIDUE_ROLE = [
 ("Tyr70",  "Y72",  "Y", "peripheral anionic site, gorge rim"),
 ("Asp72",  "D74",  "D", "gorge entrance; electrostatic steering of the cationic substrate"),
 ("Trp84",  "W86",  "W", "choline-binding (anionic) subsite; cation-pi anchor and the mid-gorge basin"),
 ("Gly118", "G120", "G", "oxyanion hole"),
 ("Gly119", "G121", "G", "oxyanion hole"),
 ("Tyr121", "Y124", "Y", "gorge lining, peripheral site"),
 ("Glu199", "E202", "E", "adjacent to the catalytic serine; anionic subsite"),
 ("Ser200", "S203", "S", "CATALYTIC NUCLEOPHILE; the atom an organophosphate phosphorylates"),
 ("Ala201", "A204", "A", "oxyanion hole"),
 ("Trp279", "W286", "W", "principal peripheral anionic site residue; first contact for an incoming ligand"),
 ("Phe288", "F295", "F", "ACYL POCKET; sets substrate size selectivity. Small residues here widen the pocket (BChE-like)"),
 ("Phe290", "F297", "F", "ACYL POCKET; with Phe288 defines the steric limit on the acyl group"),
 ("Glu327", "E334", "E", "catalytic triad (acid)"),
 ("Phe330", "F338", "F", "mid-gorge lining; cation-pi partner during descent"),
 ("Tyr334", "Y341", "Y", "gorge lining, peripheral site"),
 ("His440", "H447", "H", "catalytic triad (base); Ne2 must face Ser Og"),
]


def embed_figures(root, cap_mb=1.2):
    """Base64-embed existing figures so the dashboard stays a single file."""
    import base64
    want = [
        ("BIPHASIC/fig1_gorge_landscape", "Docking energy profile along the gorge",
         "Best-decile affinity binned by depth. A minimum between the two shaded "
         "zones is the mid-gorge basin: the thermodynamic minimum of the "
         "non-covalent complex, and the reason an unconstrained search does not "
         "return the reactive pose."),
        ("BIPHASIC/fig2_two_site", "Peripheral versus acylation site affinity",
         "Per species. Points below the diagonal favour the gorge mouth."),
        ("BIPHASIC/fig3_near_attack", "Near-attack geometry at the acylation site",
         "Restricted to poses that reached the acylation site. The shaded box is "
         "the reactive arrangement; poses outside it cannot phosphorylate."),
        ("FIGURES/fig5_ssd", "Species sensitivity distribution", ""),
        ("FIGURES/fig_composite", "Phylogeny, conservation and sensitivity", ""),
        ("FIGURES/fig4_lambda", "Phylogenetic signal likelihood", ""),
        ("MODEL_FINAL/fig_forest", "Effect sizes with confidence intervals",
         "Every predictor tested. Bars crossing zero are not distinguishable "
         "from no effect."),
        ("FIGURES/fig2_active_site", "Active-site residues across species", ""),
    ]
    out = []
    for stem, title, note in want:
        for sub in ("", "folding/", "../"):
            hit = None
            for ext, mime in ((".svg", "image/svg+xml"), (".png", "image/png")):
                fp = os.path.join(root, sub + stem + ext)
                if os.path.exists(fp) and os.path.getsize(fp) < cap_mb * 1e6:
                    if hit is None or os.path.getsize(fp) < hit[2]:
                        hit = (fp, mime, os.path.getsize(fp))
            if hit:
                with open(hit[0], "rb") as fh:
                    b = base64.b64encode(fh.read()).decode()
                out.append({"title": title, "note": note,
                            "src": f"data:{hit[1]};base64,{b}",
                            "file": os.path.basename(hit[0]),
                            "kb": round(hit[2] / 1024)})
                break
    return out


# ------------------------------------------------------------------ io ------

def rd(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def find(root, *cands):
    for c in cands:
        p = os.path.join(root, c)
        if os.path.exists(p):
            return p
    return None


def tag_of(r):
    sid = r.get("Sequence_id") or r.get("sequence_id") or ""
    if "__" in sid:
        return sid.split("__")[0]
    for k in ("species", "Species_representative", "tag"):
        if r.get(k):
            return r[k].strip().replace(" ", "_")
    return ""


# ---------------------------------------------------------------- newick ----

class Node:
    __slots__ = ("name", "length", "children", "parent")
    def __init__(self):
        self.name = ""; self.length = 0.0; self.children = []; self.parent = None


def parse_newick(text):
    text = text.strip().rstrip(";"); i = 0
    def node():
        nonlocal i
        n = Node()
        if text[i] == "(":
            i += 1
            while True:
                c = node(); c.parent = n; n.children.append(c)
                if text[i] == ",": i += 1; continue
                if text[i] == ")": i += 1; break
                raise ValueError("newick parse error")
        j = i
        while i < len(text) and text[i] not in "(),:;": i += 1
        lab = text[j:i].strip().strip("'\"")
        if not n.children: n.name = lab
        if i < len(text) and text[i] == ":":
            i += 1; j = i
            while i < len(text) and (text[i].isdigit() or text[i] in ".eE+-"): i += 1
            try: n.length = float(text[j:i])
            except ValueError: n.length = 0.0
        return n
    return node()


def leaves(n):
    return [n] if not n.children else [x for c in n.children for x in leaves(c)]


def vcv(root, order):
    depth, paths = {}, {}
    def rec(n, d, path):
        depth[n] = d
        if not n.children: paths[n.name] = (path + [n], d)
        for c in n.children: rec(c, d + c.length, path + [n])
    rec(root, 0.0, [])
    k = len(order); V = np.zeros((k, k))
    for i, A in enumerate(order):
        pa, da = paths[A]
        for j, B in enumerate(order):
            pb, _ = paths[B]
            common = [x for x in pa if x in pb]
            V[i, j] = depth[common[-1]] if common else 0.0
        V[i, i] = da
    return V


def phylo_impute(V, order, y_by, lam):
    """Conditional expectation and SD for tips without an endpoint.

    Standard BM kriging: split the lambda-transformed covariance into observed
    and unobserved blocks, then condition. The interval widens with distance
    from any measured relative, which is precisely the behaviour a read-across
    needs: a species deep inside a well-sampled clade gets a tight interval, an
    isolated one gets an honest wide one."""
    obs = [i for i, t in enumerate(order) if t in y_by]
    unk = [i for i, t in enumerate(order) if t not in y_by]
    if len(obs) < 4 or not unk:
        return {}
    d = np.diag(np.diag(V)); C = d + lam * (V - d)
    Voo = C[np.ix_(obs, obs)] + 1e-9 * np.eye(len(obs))
    Vuo = C[np.ix_(unk, obs)]
    Vuu = C[np.ix_(unk, unk)]
    yo = np.array([y_by[order[i]] for i in obs])
    Vi = np.linalg.inv(Voo)
    one = np.ones((len(obs), 1))
    # .item(), not float(): numpy 2 refuses float() on a (1,1) array
    mu = ((one.T @ Vi @ yo.reshape(-1, 1)) / (one.T @ Vi @ one)).item()
    r = yo - mu
    s2 = np.asarray(r @ Vi @ r).item() / max(1, len(obs) - 1)
    mean = np.asarray(mu + Vuo @ Vi @ r).ravel()
    var = np.diag(Vuu - Vuo @ Vi @ Vuo.T) * s2
    out = {}
    for k, i in enumerate(unk):
        sd = math.sqrt(max(0.0, float(var[k])))
        out[order[i]] = {"pred": round(float(mean[k]), 3), "sd": round(sd, 3),
                         "lo": round(float(mean[k]) - 1.96 * sd, 3),
                         "hi": round(float(mean[k]) + 1.96 * sd, 3)}
    return out


# ------------------------------------------------------------------ main ----

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tree", default=None)
    ap.add_argument("--title", default="MIE species susceptibility screening")
    ap.add_argument("--chemical", default="chlorpyrifos")
    ap.add_argument("--target", default="acetylcholinesterase")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    R = a.root

    P = {
        "endpoints": find(R, "REPORT/analysis_set.tsv", "../REPORT/analysis_set.tsv"),
        "attrition": find(R, "REPORT/table1_attrition.tsv", "../REPORT/table1_attrition.tsv"),
        "qc": find(R, "QC/model_qc.tsv", "folding/QC/model_qc.tsv"),
        "prep": find(R, "RECEPTORS/receptor_prep.tsv", "folding/RECEPTORS/receptor_prep.tsv"),
        "desc": find(R, "STRUCT/descriptors_species.tsv", "folding/STRUCT/descriptors_species.tsv"),
        "desc_icc": find(R, "STRUCT/variance_partition.tsv", "folding/STRUCT/variance_partition.tsv"),
        "biphasic": find(R, "BIPHASIC/biphasic_species.tsv", "folding/BIPHASIC/biphasic_species.tsv"),
        "pen": find(R, "BIPHASIC/penetration_descriptors.tsv", "folding/BIPHASIC/penetration_descriptors.tsv"),
        "pen_icc": find(R, "BIPHASIC/penetration_icc.tsv", "folding/BIPHASIC/penetration_icc.tsv"),
        "univ": find(R, "MODEL_FINAL/univariate_results.tsv", "MODEL/univariate_results.tsv",
                     "folding/MODEL_FINAL/univariate_results.tsv"),
        "phylo": find(R, "MODEL_FINAL/phylo_signal.tsv", "MODEL/phylo_signal.tsv",
                      "folding/MODEL_FINAL/phylo_signal.tsv"),
        "poses": find(R, "BIPHASIC/poses_all.tsv", "folding/BIPHASIC/poses_all.tsv"),
        "affvar": find(R, "DOCKANALYSIS/affinity_variance.tsv", "folding/DOCKANALYSIS/affinity_variance.tsv"),
    }
    # A DATED SPECIES TREE IS STRONGLY PREFERRED. Its branch lengths are
    # divergence times; a gene tree's are AChE sequence divergence, so
    # imputation on a gene tree borrows strength according to how fast the
    # protein evolved rather than how long ago the species split.
    tree_p = a.tree or find(R,
        "TIMETREE/species_tree_reconciled.nwk",
        "folding/TIMETREE/species_tree_reconciled.nwk",
        "../TIMETREE/species_tree_reconciled.nwk",
        "PHYLO/ache.raxml.support", "folding/PHYLO/ache.raxml.support",
        "TREE/ache.tree", "folding/TREE/ache.tree")
    print("inputs found:")
    for k, v in P.items():
        print(f"  {k:10s} {v or 'MISSING'}")
    print(f"  {'tree':10s} {tree_p or 'MISSING'}")

    # ---- species records --------------------------------------------------
    S = defaultdict(lambda: {"tiers": {}})
    for r in rd(P["endpoints"]):
        t = tag_of(r)
        if not t: continue
        S[t].update(tag=t,
                    name=r.get("Species_representative") or t.replace("_", " "),
                    klass=r.get("Class", ""), order=r.get("Order", ""),
                    family=r.get("Family", ""),
                    lc50=fnum(r.get("log10_LC50")),
                    members=r.get("member_species", ""),
                    n_collapsed=r.get("n_species_collapsed", ""),
                    seq_source=r.get("Sequence_source", ""),
                    clade=r.get("AChE_clade", ""), tier_flag=r.get("Tier", ""))
    for r in rd(P["qc"]):
        t = r.get("species", "")
        if not t: continue
        S[t].setdefault("tag", t); S[t].setdefault("name", t.replace("_", " "))
        S[t].update(plddt=fnum(r.get("plddt_site_mean")),
                    plddt_crit=fnum(r.get("plddt_critical_min")),
                    seed_rmsd=fnum(r.get("seed_rmsd_site_allatom_mean")),
                    qc=r.get("verdict", ""),
                    active_site="".join(r.get(k, "-") for k in
                                        ("Tyr70", "Asp72", "Trp84", "Gly118",
                                         "Gly119", "Tyr121", "Glu199", "Ser200",
                                         "Ala201", "Trp279", "Phe288", "Phe290",
                                         "Glu327", "Phe330", "Tyr334", "His440")))
    for r in rd(P["prep"]):
        t = r.get("species", "")
        if t in S:
            S[t].update(gorge=fnum(r.get("gorge_length")),
                        triad=r.get("triad_ok", ""),
                        his_flip=r.get("his_flipped", ""))
    DESC_KEYS = []
    for r in rd(P["desc"]):
        t = r.get("species", "")
        if not t: continue
        d = {}
        for k, v in r.items():
            if k in ("species", "n_models") or k.endswith("_sd"): continue
            fv = fnum(v)
            if fv is not None:
                d[k] = fv
                if k not in DESC_KEYS: DESC_KEYS.append(k)
        S[t].setdefault("tag", t); S[t].setdefault("name", t.replace("_", " "))
        S[t]["desc"] = d
    for r in rd(P["pen"]):
        t = r.get("species", "")
        if not t: continue
        p = {k: fnum(v) for k, v in r.items()
             if k.startswith("pen_") and not k.endswith("_sd") and fnum(v) is not None}
        S[t].setdefault("tag", t); S[t].setdefault("name", t.replace("_", " "))
        S[t]["pen"] = p
    for r in rd(P["biphasic"]):
        t = r.get("species", "")
        if t in S and r.get("ligand") == "chlorpyrifos_oxon":
            S[t].update(dG_A=fnum(r.get("dG_A")), dG_M=fnum(r.get("dG_M")),
                        dG_P=fnum(r.get("dG_P")), ddG=fnum(r.get("ddG_PA")))

    # ---- model status: this gates the prediction tab ----------------------
    univ = rd(P["univ"])
    phy = rd(P["phylo"])
    lam = fnum(phy[0].get("lambda")) if phy else None
    lam_lr = fnum(phy[0].get("LR")) if phy else None
    surv, surv_raw = [], []
    for r in univ:
        q = fnum(r.get("pgls_fix_p_holm"))
        if q is not None and q < 0.05:
            surv.append(r["predictor"])
        q2 = fnum(r.get("p_holm"))
        if q2 is not None and q2 < 0.05:
            surv_raw.append(r["predictor"])

    # ---- phylogenetic prediction -----------------------------------------
    pred = {}
    tree_note = "no tree available"
    if tree_p and os.path.exists(tree_p):
        try:
            root = parse_newick(open(tree_p).read())
            lf = leaves(root)
            gene_like = sum(1 for l in lf if "__orthodb__" in l.name)
            m, seen = {}, set()
            for l in lf:
                for cand in (l.name, l.name.split("__")[0]):
                    if cand in S:
                        # A gene tree carries several sequences per species and
                        # they all rename to the same tag. Left unchecked the
                        # tip order contains duplicates and the covariance
                        # matrix is singular. Keep the first tip per species.
                        if cand in seen:
                            break
                        seen.add(cand); m[l.name] = cand
                        break
            for l in lf:
                l.name = m.get(l.name, "__drop__" + l.name)
            order = [l.name for l in leaves(root) if l.name in S]
            n_dupe = sum(1 for l in lf if l.name.startswith("__drop__")
                         and l.name[8:].split("__")[0] in S)
            if len(order) >= 6:
                V = vcv(root, order)
                y_by = {t: S[t]["lc50"] for t in order
                        if S[t].get("lc50") is not None}
                pred = phylo_impute(V, order, y_by, lam if lam is not None else 1.0)
                kind = ("GENE tree (branch lengths are AChE sequence "
                        "divergence, not time)" if gene_like > 0
                        else "dated species tree")
                tree_note = (f"{os.path.basename(tree_p)} \u2014 {kind}; "
                             f"{len(order)} species, {len(y_by)} measured")
                if gene_like > 0:
                    print("\n  *** WARNING: predictions are on a GENE tree ***")
                    print("  Branch lengths are AChE sequence divergence, not "
                          "divergence time, so a fast-evolving lineage looks "
                          "further from its relatives than it is and its "
                          "prediction interval is inflated accordingly.")
                    print("  Submit ALL species to TimeTree, rerun "
                          "timetree_helper.py, and pass --tree "
                          "TIMETREE/species_tree_reconciled.nwk.")
                if n_dupe:
                    print(f"  ({n_dupe} additional tips mapped to species "
                          "already present and were dropped; one sequence per "
                          "species is used)")
        except Exception as e:
            tree_note = f"tree could not be used: {e}"
    for t, v in pred.items():
        S[t]["pred"] = v
    print(f"\nphylogenetic prediction: {len(pred)} species imputed ({tree_note})")

    # ---- pose zone occupancy ---------------------------------------------
    zones = {}
    for r in rd(P["poses"]):
        lig, z = r.get("ligand", ""), r.get("zone", "")
        if lig:
            zones.setdefault(lig, Counter())[z] += 1
    zone_tbl = []
    for lig, c in sorted(zones.items()):
        n = sum(c.values())
        zone_tbl.append({"ligand": lig, "n": n,
                         "A": round(100 * c["A"] / n, 1),
                         "M": round(100 * c["M"] / n, 1),
                         "P": round(100 * c["P"] / n, 1)})

    # per-species pose-depth histogram: the shape of where each ligand sits
    NB = 12
    hist = defaultdict(lambda: defaultdict(lambda: [0] * NB))
    for r in rd(P["poses"]):
        t, lig = r.get("species", ""), r.get("ligand", "")
        f = fnum(r.get("frac_depth"))
        if not t or not lig or f is None:
            continue
        b = min(NB - 1, max(0, int(f / 1.1 * NB)))
        hist[t][lig][b] += 1
    for t, d in hist.items():
        if t in S:
            S[t]["depth"] = {k: v for k, v in d.items()}

    figs = embed_figures(R)
    print(f"\nfigures embedded: {len(figs)}"
          + ("" if figs else " (none found; run the figure-producing scripts)"))
    for f in figs:
        print(f"    {f['file']:32s} {f['kb']:5d} KB")

    icc_tbl = [{"descriptor": r.get("descriptor") or r.get("ligand", ""),
                "ICC": fnum(r.get("ICC")), "usable": r.get("usable", "")}
               for r in rd(P["desc_icc"]) + rd(P["pen_icc"]) + rd(P["affvar"])]
    icc_tbl = [r for r in icc_tbl if r["ICC"] is not None]

    payload = {
        "meta": {"title": a.title, "chemical": a.chemical, "target": a.target,
                 "generated": time.strftime("%Y-%m-%d %H:%M"),
                 "tree_note": tree_note},
        "status": {"lambda": lam, "lambda_LR": lam_lr,
                   "n_tested": len(univ), "surviving": surv,
                   "surviving_raw": surv_raw,
                   "n_species_endpoint": sum(1 for v in S.values()
                                             if v.get("lc50") is not None),
                   "n_species_total": len(S),
                   "n_predicted": len(pred)},
        "attrition": rd(P["attrition"]),
        "species": sorted(S.values(), key=lambda v: (v.get("klass") or "zz",
                                                     v.get("name") or "")),
        "descriptor_keys": DESC_KEYS,
        "univariate": univ, "icc": icc_tbl, "zones": zone_tbl,
        "figures": figs, "residues": RESIDUE_ROLE,
    }

    out = os.path.join(a.outdir, "dashboard.html")
    with open(out, "w") as fh:
        fh.write(HTML.replace("__DATA__", json.dumps(payload))
                     .replace("__TITLE__", html.escape(a.title)))
    with open(os.path.join(a.outdir, "dashboard_data.json"), "w") as fh:
        json.dump(payload, fh, indent=1)

    print(f"\n  {out}")
    print(f"  {a.outdir}/dashboard_data.json")
    print(f"\n  species: {payload['status']['n_species_total']} "
          f"({payload['status']['n_species_endpoint']} measured, "
          f"{len(pred)} phylogenetically predicted)")
    if surv:
        print(f"  predictors surviving validation: {', '.join(surv)}")
    else:
        print("  predictors surviving validation: NONE. The dashboard states "
              "this in its banner and does not offer a structure-based "
              "ranking; prediction falls back to phylogenetic read-across, "
              "which the measured lambda supports.")


HTML = r"""<!doctype html><html lang="en"><meta charset="utf-8">
<title>__TITLE__</title>
<style>
:root{--bg:#f6f7f9;--card:#fff;--ink:#1c2530;--mut:#63707e;--line:#dde3ea;
--ok:#1a7f4b;--warn:#b8770c;--bad:#b3261e;--acc:#2f4b7c}
*{box-sizing:border-box}body{margin:0;font:14px/1.5 -apple-system,Segoe UI,Roboto,Helvetica,sans-serif;
background:var(--bg);color:var(--ink)}
header{background:var(--acc);color:#fff;padding:16px 22px}
header h1{margin:0;font-size:19px;font-weight:600}
header .sub{opacity:.85;font-size:12.5px;margin-top:3px}
.banner{padding:12px 22px;font-size:13.5px;border-bottom:1px solid var(--line)}
.banner.bad{background:#fdecea;color:#7a1b14;border-left:5px solid var(--bad)}
.banner.ok{background:#e9f6ef;color:#125c37;border-left:5px solid var(--ok)}
nav{display:flex;gap:2px;background:#eceff3;padding:0 14px;border-bottom:1px solid var(--line);
flex-wrap:wrap}
nav button{border:0;background:none;padding:11px 15px;cursor:pointer;font-size:13.5px;
color:var(--mut);border-bottom:3px solid transparent}
nav button.on{color:var(--acc);border-bottom-color:var(--acc);font-weight:600}
main{padding:18px 22px;max-width:1500px}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:16px;
margin-bottom:16px}
.card h2{margin:0 0 4px;font-size:15px}
.card p.note{margin:0 0 12px;color:var(--mut);font-size:12.5px}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{padding:6px 9px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}
th{background:#f2f5f8;cursor:pointer;position:sticky;top:0;font-weight:600}
tbody tr:hover{background:#f4f8fd;cursor:pointer}
.wrap{max-height:560px;overflow:auto;border:1px solid var(--line);border-radius:6px}
.pill{padding:1px 8px;border-radius:10px;font-size:11px;font-weight:600}
.p-ok{background:#e2f3ea;color:var(--ok)}.p-warn{background:#fdf1dc;color:var(--warn)}
.p-bad{background:#fdecea;color:var(--bad)}.p-mut{background:#eef1f4;color:var(--mut)}
input,select{padding:6px 9px;border:1px solid var(--line);border-radius:5px;font-size:13px}
.kv{display:grid;grid-template-columns:190px 1fr;gap:4px 12px;font-size:12.5px}
.kv div:nth-child(odd){color:var(--mut)}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:12px}
.stat{background:#f2f5f8;border-radius:6px;padding:11px}
.stat b{display:block;font-size:21px;color:var(--acc)}
.stat span{font-size:11.5px;color:var(--mut)}
#detail{position:fixed;right:0;top:0;height:100%;width:430px;background:#fff;
border-left:1px solid var(--line);box-shadow:-3px 0 14px rgba(0,0,0,.08);padding:18px;
overflow:auto;transform:translateX(100%);transition:.18s}
#detail.on{transform:none}
#detail .x{float:right;cursor:pointer;color:var(--mut);font-size:19px}
code{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:11.5px;
background:#f2f5f8;padding:1px 4px;border-radius:3px}
.mono{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:11.5px}
</style>
<header>
  <h1 id="ttl"></h1><div class="sub" id="sub"></div>
</header>
<div id="banner"></div>
<nav id="nav"></nav>
<main id="main"></main>
<div id="detail"></div>
<script>
const D = __DATA__;
const $ = s => document.querySelector(s);
const esc = s => String(s==null?"":s).replace(/[&<>]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;"}[c]));
const fx = (v,n=2) => (v==null||v===""||isNaN(v))?"":Number(v).toFixed(n);

$("#ttl").textContent = D.meta.title;
$("#sub").textContent = `${D.meta.chemical} / ${D.meta.target} · generated ${D.meta.generated}`;

// ---- banner: the validation gate ----
(function(){
  const s=D.status, ok=s.surviving.length>0;
  const b=$("#banner"); b.className="banner "+(ok?"ok":"bad");
  b.innerHTML = ok
   ? `<b>${s.surviving.length} predictor(s) survived phylogenetic correction and multiple testing:</b> ${s.surviving.map(esc).join(", ")}. Structure-based ranking is supported for these predictors only.`
   : `<b>No predictor survived validation.</b> ${s.n_tested} predictors were tested; ${s.surviving_raw.length} were significant before correction, none after phylogenetic correction and multiple testing. Phylogenetic signal in the endpoint is &lambda; = ${fx(s.lambda,3)}. <u>This dashboard therefore does not offer a structure-based susceptibility ranking.</u> Prediction below is phylogenetic read-across, which the measured &lambda; supports.`;
})();

const TABS=[
 ["t0","Tier 0 · Coverage",tier0],["t1","Tier 1 · Target",tier1],
 ["t2","Tier 2 · Structure",tier2],["t3","Tier 3 · Binding",tier3],
 ["t4","Tier 4 · Inference",tier4],["t5","Tier 5 · Prediction",tier5],
 ["mx","Mechanism",mech],["fg","Figures",figs],
 ["sp","Species browser",browser]];
TABS.forEach(([id,lab])=>{const b=document.createElement("button");b.textContent=lab;
 b.onclick=()=>show(id);b.id="b_"+id;$("#nav").appendChild(b);});
function show(id){TABS.forEach(([i])=>$("#b_"+i).classList.toggle("on",i===id));
 $("#main").innerHTML="";TABS.find(t=>t[0]===id)[2]($("#main"));}

function card(p,t,note,inner){const d=document.createElement("div");d.className="card";
 d.innerHTML=`<h2>${t}</h2>${note?`<p class="note">${note}</p>`:""}${inner||""}`;
 p.appendChild(d);return d;}

function table(rows,cols,onclick){
 if(!rows.length) return "<p class='note'>no data</p>";
 let h="<div class='wrap'><table><thead><tr>"+cols.map(c=>`<th>${esc(c.t||c.k)}</th>`).join("")+"</tr></thead><tbody>";
 rows.forEach((r,i)=>{h+=`<tr data-i='${i}'>`+cols.map(c=>`<td>${c.f?c.f(r):esc(r[c.k])}</td>`).join("")+"</tr>";});
 return h+"</tbody></table></div>";}

function wire(el,rows,cb){
 el.querySelectorAll("tbody tr").forEach(tr=>tr.onclick=()=>cb(rows[+tr.dataset.i]));
 el.querySelectorAll("th").forEach((th,ci)=>th.onclick=()=>{
  const tb=th.closest("table").tBodies[0];
  const rs=[...tb.rows]; const dir=th.dataset.d==="1"?-1:1; th.dataset.d=dir===1?"1":"";
  rs.sort((a,b)=>{const x=a.cells[ci].textContent,y=b.cells[ci].textContent;
   const nx=parseFloat(x),ny=parseFloat(y);
   return (!isNaN(nx)&&!isNaN(ny))?(nx-ny)*dir:x.localeCompare(y)*dir;});
  rs.forEach(r=>tb.appendChild(r));});}

function pill(txt,cls){return `<span class="pill ${cls}">${esc(txt)}</span>`;}
function iccPill(v){v=+v;return pill(fx(v,2), v>=0.75?"p-ok":(v>=0.5?"p-warn":"p-bad"));}

// ---- tiers ----
function tier0(p){
 const s=D.status;
 card(p,"Assessment coverage","How many species could be assessed at all, and where the rest were lost. This is the first thing a regulator needs: the tool's applicability domain.",
  `<div class="grid">
   <div class="stat"><b>${s.n_species_total}</b><span>species with a verified target</span></div>
   <div class="stat"><b>${s.n_species_endpoint}</b><span>with a measured endpoint</span></div>
   <div class="stat"><b>${s.n_predicted}</b><span>predicted by read-across</span></div>
   <div class="stat"><b>${fx(s.lambda,3)}</b><span>Pagel's &lambda; in the endpoint</span></div></div>`);
 if(D.attrition.length){
  const c=card(p,"Attrition","Every species lost, and why. A gap here is a data gap in the world, not a failure of the tool.","");
  const cols=Object.keys(D.attrition[0]).map(k=>({k}));
  c.insertAdjacentHTML("beforeend",table(D.attrition,cols));}
}
function tier1(p){
 const rows=D.species.filter(x=>x.clade||x.qc);
 const c=card(p,"Target verification","Each species' sequence was assigned to an AChE clade by phylogenetic placement, not by annotation keywords. Tier B/C reflects congener proxies or multiple candidates.","");
 c.insertAdjacentHTML("beforeend",table(rows,[
  {k:"name",t:"species"},{k:"klass",t:"class"},{k:"order"},{k:"family"},
  {k:"clade",t:"AChE clade"},{k:"seq_source",t:"sequence source"},
  {k:"tier_flag",t:"tier",f:r=>r.tier_flag?pill(r.tier_flag,r.tier_flag==="A"?"p-ok":(r.tier_flag==="B"?"p-warn":"p-bad")):""},
  {k:"n_collapsed",t:"species collapsed"}]));
 wire(c,rows,detail);
}
function tier2(p){
 const rows=D.species.filter(x=>x.plddt!=null);
 const c=card(p,"Structure quality and geometry","Model confidence is reported AT THE ACTIVE SITE, not as a whole-model mean: a model can average 92 with a ragged acyl pocket. Seed RMSD is the noise floor from five independent predictions.","");
 c.insertAdjacentHTML("beforeend",table(rows,[
  {k:"name",t:"species"},{k:"klass",t:"class"},
  {k:"plddt",t:"pLDDT active site",f:r=>fx(r.plddt,1)},
  {k:"plddt_crit",t:"pLDDT critical min",f:r=>{const v=+r.plddt_crit;
    return pill(fx(v,1), v>=90?"p-ok":(v>=70?"p-warn":"p-bad"));}},
  {k:"seed_rmsd",t:"seed RMSD (Å)",f:r=>{const v=+r.seed_rmsd;
    return pill(fx(v,2), v<=0.5?"p-ok":(v<=1.0?"p-warn":"p-bad"));}},
  {k:"gorge",t:"gorge length (Å)",f:r=>fx(r.gorge,2)},
  {k:"triad",t:"triad ok"},
  {k:"active_site",t:"active site",f:r=>`<span class="mono">${esc(r.active_site||"")}</span>`}]));
 wire(c,rows,detail);
 if(D.icc.length){
  const c2=card(p,"Descriptor reliability (ICC)","Fraction of each descriptor's variance that is between species rather than between replicate predictions. Below 0.5 the descriptor is mostly noise and its regression coefficient is not interpretable.","");
  c2.insertAdjacentHTML("beforeend",table(D.icc,[
   {k:"descriptor"},{k:"ICC",f:r=>iccPill(r.ICC)},{k:"usable"}]));}
}
function tier3(p){
 if(D.zones.length){
  card(p,"Where ligands sit in the gorge","Percentage of all docked poses in each zone: A = acylation site at the base, M = mid-gorge basin, P = peripheral site at the mouth. A large ligand trapped mid-gorge while the natural substrate reaches the base is the two-site mechanism, recovered without being assumed.",
   table(D.zones,[{k:"ligand"},{k:"n",t:"poses"},
    {k:"A",t:"% acylation site"},{k:"M",t:"% mid-gorge"},{k:"P",t:"% peripheral"}]));}
 const rows=D.species.filter(x=>x.dG_A!=null||x.dG_M!=null||x.dG_P!=null);
 if(rows.length){
  const c=card(p,"Two-site affinities","Best affinity per zone. These are non-covalent scores: organophosphate potency is k2/Kd and nothing here addresses k2. Read the contrasts, not the absolute values.","");
  c.insertAdjacentHTML("beforeend",table(rows,[
   {k:"name",t:"species"},{k:"klass",t:"class"},
   {k:"dG_A",t:"ΔG acylation",f:r=>fx(r.dG_A)},{k:"dG_M",t:"ΔG mid-gorge",f:r=>fx(r.dG_M)},
   {k:"dG_P",t:"ΔG peripheral",f:r=>fx(r.dG_P)},{k:"ddG",t:"ΔΔG P−A",f:r=>fx(r.ddG)}]));
  wire(c,rows,detail);}
}
function tier4(p){
 const u=D.univariate;
 card(p,"Phylogenetic signal in the endpoint",
  "Pagel's &lambda; for the endpoint with no predictor. Near 1 means sensitivity is largely inherited, so a raw correlation between any taxon-linked trait and sensitivity is expected whether or not it is causal.",
  `<div class="grid"><div class="stat"><b>${fx(D.status.lambda,3)}</b><span>&lambda;</span></div>
   <div class="stat"><b>${fx(D.status.lambda_LR,2)}</b><span>LR vs &lambda;=0 (>3.84 rejects)</span></div>
   <div class="stat"><b>${D.status.surviving_raw.length}/${D.status.n_tested}</b><span>significant before correction</span></div>
   <div class="stat"><b>${D.status.surviving.length}/${D.status.n_tested}</b><span>significant after correction</span></div></div>
   <p class="note" style="margin-top:10px">Tree: ${esc(D.meta.tree_note)}</p>`);
 if(u.length){
  const c=card(p,"Every predictor tested","Reported whether or not significant. <code>p_holm</code> is uncorrected for phylogeny; <code>pgls_fix_p_holm</code> is the corrected test at fixed &lambda;, which is the one to read. A &lambda; at the 0 or 1 boundary means the correction degenerated and that row's estimated-&lambda; p is not a correction.","");
  c.insertAdjacentHTML("beforeend",table(u,[
   {k:"predictor"},{k:"family"},{k:"n"},
   {k:"pearson_r",t:"r",f:r=>fx(r.pearson_r,3)},
   {k:"ICC",f:r=>r.ICC?iccPill(r.ICC):""},
   {k:"p_holm",t:"p (raw, Holm)",f:r=>{const v=+r.p_holm;
     return isNaN(v)?"":pill(v.toFixed(3),v<0.05?"p-warn":"p-mut");}},
   {k:"pgls_lambda",t:"λ resid"},
   {k:"pgls_fix_p_holm",t:"p (corrected, Holm)",f:r=>{const v=+r.pgls_fix_p_holm;
     return isNaN(v)?"":pill(v.toFixed(3),v<0.05?"p-ok":"p-mut");}}]));}
}
function tier5(p){
 const gate=D.status.surviving.length>0;
 card(p,"How these predictions are made",
  gate?"Structure-based predictors survived validation and may be used alongside read-across."
      :"No structural or docking predictor survived validation, so susceptibility is <b>not</b> predicted from protein geometry here. Instead it is imputed phylogenetically: the conditional expectation of an unmeasured species given its relatives and the tree, under Brownian motion with the measured &lambda;. The interval widens with distance from any measured relative, so an isolated species is reported with an honestly wide interval rather than a confident guess.","");
 const rows=D.species.filter(x=>x.pred).map(x=>Object.assign({},x,
   {p_pred:x.pred.pred,p_lo:x.pred.lo,p_hi:x.pred.hi,p_sd:x.pred.sd,
    width:+(x.pred.hi-x.pred.lo).toFixed(2)}));
 if(!rows.length){card(p,"No predictions","No tree was available, or every species already has a measured endpoint.","");return;}
 const c=card(p,`Predicted susceptibility (${rows.length} species)`,
  "log10 LC50 in mg/L. Narrower intervals mean the species sits inside a well-sampled clade. Sort by interval width to see which predictions are worth acting on.","");
 c.insertAdjacentHTML("beforeend",table(rows,[
  {k:"name",t:"species"},{k:"klass",t:"class"},{k:"order"},
  {k:"p_pred",t:"predicted log10 LC50",f:r=>fx(r.p_pred,2)},
  {k:"p_lo",t:"95% lower",f:r=>fx(r.p_lo,2)},
  {k:"p_hi",t:"95% upper",f:r=>fx(r.p_hi,2)},
  {k:"width",t:"interval width",f:r=>{const v=+r.width;
    return pill(v.toFixed(2), v<2?"p-ok":(v<4?"p-warn":"p-bad"));}},
  {k:"pred_mgL",t:"LC50 mg/L",f:r=>Math.pow(10,r.p_pred).toPrecision(3)}]));
 wire(c,rows,detail);
}
function browser(p){
 const c=card(p,"All species","Click any row for the full record: taxonomy, target assignment, model quality, active-site residues, descriptors, binding and prediction.",
  `<div style="margin-bottom:10px"><input id="q" placeholder="filter by name, class, order…" style="width:320px"></div><div id="tw"></div>`);
 const draw=f=>{const rows=D.species.filter(x=>!f||
   JSON.stringify([x.name,x.klass,x.order,x.family]).toLowerCase().includes(f.toLowerCase()));
  $("#tw").innerHTML=table(rows,[
   {k:"name",t:"species"},{k:"klass",t:"class"},{k:"order"},{k:"family"},
   {k:"lc50",t:"measured log10 LC50",f:r=>r.lc50!=null?fx(r.lc50,2):pill("predicted","p-mut")},
   {k:"qc",t:"model QC",f:r=>r.qc?pill(r.qc,r.qc==="pass"?"p-ok":"p-bad"):""},
   {k:"clade",t:"clade"}]);
  wire($("#tw"),rows,detail);};
 draw(""); $("#q").oninput=e=>draw(e.target.value);
}
function gorgeSVG(){
 const z=D.zones||[]; const g=l=>z.find(x=>x.ligand===l)||{A:0,M:0,P:0};
 const ox=g("chlorpyrifos_oxon"), ach=g("acetylcholine"), tcp=g("TCP");
 const bar=(v,x,y,c)=>`<rect x="${x}" y="${y}" width="${Math.max(1,v*1.15)}" height="11" fill="${c}" rx="2"/>
   <text x="${x+v*1.15+5}" y="${y+9}" font-size="10" fill="#63707e">${v}%</text>`;
 return `<svg viewBox="0 0 720 330" style="width:100%;max-width:720px">
 <defs><linearGradient id="gg" x1="0" y1="0" x2="0" y2="1">
   <stop offset="0" stop-color="#dbe6f2"/><stop offset="1" stop-color="#f4f7fa"/></linearGradient></defs>
 <path d="M60,28 L200,28 L152,175 L142,300 L118,300 L108,175 Z" fill="url(#gg)" stroke="#b9c6d4"/>
 <text x="130" y="20" font-size="11" text-anchor="middle" fill="#63707e">gorge mouth</text>
 <circle cx="130" cy="45" r="5" fill="#2f4b7c"/><text x="212" y="42" font-size="11.5"><tspan font-weight="600">P-site</tspan> · Trp279, Tyr70, Asp72, Tyr121</text>
 <text x="212" y="56" font-size="10.5" fill="#63707e">capture and orientation of the incoming ligand</text>
 <circle cx="130" cy="165" r="5" fill="#c0392b"/><text x="212" y="162" font-size="11.5"><tspan font-weight="600">mid-gorge basin</tspan> · Trp84, Phe330</text>
 <text x="212" y="176" font-size="10.5" fill="#63707e">cation-pi minimum; where a rigid docking search settles</text>
 <circle cx="130" cy="288" r="5" fill="#1a7f4b"/><text x="212" y="285" font-size="11.5"><tspan font-weight="600">A-site</tspan> · Ser200, His440, Glu327, Phe288/290</text>
 <text x="212" y="299" font-size="10.5" fill="#63707e">catalytic triad, acyl pocket, oxyanion hole; phosphorylation happens here</text>
 <text x="470" y="42" font-size="10.5" font-weight="600">poses at P-site</text>
 ${bar(ox.P,470,48,"#c0392b")}${bar(ach.P,470,64,"#1a7f4b")}${bar(tcp.P,470,80,"#2f4b7c")}
 <text x="470" y="162" font-size="10.5" font-weight="600">poses mid-gorge</text>
 ${bar(ox.M,470,168,"#c0392b")}${bar(ach.M,470,184,"#1a7f4b")}${bar(tcp.M,470,200,"#2f4b7c")}
 <text x="470" y="262" font-size="10.5" font-weight="600">poses at A-site</text>
 ${bar(ox.A,470,268,"#c0392b")}${bar(ach.A,470,284,"#1a7f4b")}${bar(tcp.A,470,300,"#2f4b7c")}
 <text x="614" y="20" font-size="10"><tspan fill="#c0392b">&#9632;</tspan> oxon <tspan fill="#1a7f4b">&#9632;</tspan> ACh <tspan fill="#2f4b7c">&#9632;</tspan> TCP</text>
 </svg>`;}

function rxnSVG(){
 const b=(x,t,s2,c)=>`<rect x="${x}" y="34" width="126" height="46" rx="6" fill="${c}" stroke="#c8d3de"/>
  <text x="${x+63}" y="55" font-size="11.5" text-anchor="middle" font-weight="600">${t}</text>
  <text x="${x+63}" y="70" font-size="10" text-anchor="middle" fill="#63707e">${s2}</text>`;
 const ar=(x,l)=>`<line x1="${x}" y1="57" x2="${x+34}" y2="57" stroke="#8b98a6" stroke-width="1.4" marker-end="url(#ah)"/>
  <text x="${x+17}" y="49" font-size="9.5" text-anchor="middle" fill="#63707e">${l}</text>`;
 return `<svg viewBox="0 0 740 110" style="width:100%;max-width:740px">
 <defs><marker id="ah" markerWidth="7" markerHeight="7" refX="6" refY="3" orient="auto">
  <path d="M0,0 L6,3 L0,6 z" fill="#8b98a6"/></marker></defs>
 ${b(4,"chlorpyrifos","P=S · not an inhibitor","#f7f9fb")}${ar(130,"CYP450")}
 ${b(164,"chlorpyrifos-oxon","P=O · the inhibitor","#fdecea")}${ar(290,"capture")}
 ${b(324,"Michaelis complex","non-covalent, reversible","#f7f9fb")}${ar(450,"k2")}
 ${b(484,"phosphorylated Ser","covalent, inactive","#e9f6ef")}${ar(610,"aging")}
 <text x="648" y="57" font-size="11" font-weight="600">aged adduct</text>
 <text x="648" y="71" font-size="10" fill="#63707e">irreversible</text>
 <text x="4" y="100" font-size="10.5" fill="#63707e">Potency is k&#8322;/K&#7480;. Docking addresses K&#7480; only; nothing in a docking score touches k&#8322;.</text>
 </svg>`;}

function mech(p){
 card(p,"Where the chemical goes, and why that matters",
  "The active site is a narrow gorge with two binding sites. Organophosphates are captured at the peripheral site, descend past a cation-pi basin, and phosphorylate the catalytic serine at the base. The bars are measured from every docked pose in this run, not assumed.",
  gorgeSVG());
 card(p,"Activation and inhibition pathway",
  "Chlorpyrifos itself cannot phosphorylate the enzyme; it requires CYP450 desulfuration to the oxon. It is included in this pipeline as a negative control precisely for that reason: if the parent predicts toxicity as well as the oxon, the pipeline is ranking lipophilicity rather than target engagement.",
  rxnSVG());
 const zc=card(p,"Gorge penetration by ligand",
  "Percentage of all poses in each zone. A bulky phosphotriester trapped mid-gorge while the small natural substrate reaches the base is size- and moiety-dependent exclusion at the constriction.","");
 zc.insertAdjacentHTML("beforeend",table(D.zones,[
  {k:"ligand"},{k:"n",t:"poses"},
  {k:"A",t:"% acylation site",f:r=>`${r.A}% ${sbar(r.A,"#1a7f4b")}`},
  {k:"M",t:"% mid-gorge",f:r=>`${r.M}% ${sbar(r.M,"#c0392b")}`},
  {k:"P",t:"% peripheral",f:r=>`${r.P}% ${sbar(r.P,"#2f4b7c")}`}]));
 const rc=card(p,"What each diagnostic residue does",
  "Torpedo mature numbering, with the mammalian equivalent. These are the positions read for every species; the reference column is the canonical residue.","");
 rc.insertAdjacentHTML("beforeend",table(D.residues.map(r=>({
   pos:r[0],mam:r[1],ref:r[2],role:r[3]})),
   [{k:"pos",t:"position"},{k:"mam",t:"mammalian"},{k:"ref",t:"canonical"},
    {k:"role",t:"function",f:r=>`<span style="white-space:normal">${esc(r.role)}</span>`}]));
}
function sbar(v,c){return `<span style="display:inline-block;height:7px;width:${Math.max(1,v*0.7)}px;background:${c};border-radius:2px;vertical-align:middle"></span>`;}

function figs(p){
 if(!D.figures.length){card(p,"No figures found","Run the figure-producing scripts, then rebuild the dashboard.","");return;}
 D.figures.forEach(f=>card(p,f.title,f.note||"",
  `<img src="${f.src}" style="width:100%;max-width:1100px;border:1px solid var(--line);border-radius:6px">
   <p class="note" style="margin-top:6px">${esc(f.file)} · ${f.kb} KB, embedded</p>`));
}

function asTable(r){
 const a=(r.active_site||"").split("");
 if(!a.length) return "";
 let h=`<p class="note" style="margin:8px 0 4px">Active site, Torpedo numbering. Orange marks a substitution from the canonical residue.</p>
  <table style="font-size:11.5px"><tbody>`;
 D.residues.forEach((rr,i)=>{const got=a[i]||"-";const dif=got!==rr[2]&&got!=="-";
  h+=`<tr><td class="mono">${esc(rr[0])}</td>
   <td class="mono" style="font-weight:600;color:${dif?"#b3261e":"#1c2530"}">${esc(got)}</td>
   <td class="mono" style="color:#63707e">${esc(rr[2])}</td>
   <td style="white-space:normal;font-size:10.5px;color:${dif?"#b3261e":"#63707e"}">${esc(rr[3])}</td></tr>`;});
 return h+"</tbody></table>";}

function spark(hist,col){
 if(!hist) return "";
 const mx=Math.max(...hist,1), w=13, h=34;
 let s="";
 hist.forEach((v,i)=>{const bh=v/mx*h;
  s+=`<rect x="${i*w}" y="${h-bh}" width="${w-2}" height="${bh}" fill="${col}" rx="1"/>`;});
 return `<svg viewBox="0 0 ${hist.length*w} ${h+13}" style="width:100%;max-width:340px">
  ${s}<text x="0" y="${h+11}" font-size="9" fill="#63707e">A-site</text>
  <text x="${hist.length*w}" y="${h+11}" font-size="9" text-anchor="end" fill="#63707e">P-site</text></svg>`;}

function detail(r){
 const d=$("#detail"); d.classList.add("on");
 const kv=o=>Object.entries(o).map(([k,v])=>`<div>${esc(k)}</div><div>${esc(v)}</div>`).join("");
 let h=`<span class="x" onclick="document.getElementById('detail').classList.remove('on')">&times;</span>
  <h2 style="margin:0 0 3px"><i>${esc(r.name)}</i></h2>
  <p class="note">${esc([r.klass,r.order,r.family].filter(Boolean).join(" · "))}</p>`;
 h+=`<h3 style="font-size:13px">Endpoint</h3><div class="kv">${kv({
  "measured log10 LC50": r.lc50!=null?fx(r.lc50,3):"none",
  "LC50 mg/L": r.lc50!=null?Math.pow(10,r.lc50).toPrecision(3):"",
  "species collapsed into this row": r.n_collapsed||"",
  "members": r.members||"", "sequence source": r.seq_source||"",
  "evidence tier": r.tier_flag||""})}</div>`;
 if(r.pred) h+=`<h3 style="font-size:13px">Phylogenetic prediction</h3><div class="kv">${kv({
  "predicted log10 LC50": r.pred.pred, "95% interval": `${r.pred.lo} to ${r.pred.hi}`,
  "SD": r.pred.sd, "LC50 mg/L": Math.pow(10,r.pred.pred).toPrecision(3)})}</div>`;
 if(r.plddt!=null) h+=`<h3 style="font-size:13px">Model</h3><div class="kv">${kv({
  "QC verdict": r.qc||"", "pLDDT active site": fx(r.plddt,1),
  "pLDDT critical minimum": fx(r.plddt_crit,1),
  "between-seed RMSD (Å)": fx(r.seed_rmsd,2),
  "catalytic triad intact": r.triad||"", "histidine flipped": r.his_flip||""})}</div>
  ${asTable(r)}`;
 if(r.desc) h+=`<h3 style="font-size:13px">Structural descriptors</h3><div class="kv">${kv(r.desc)}</div>`;
 if(r.dG_A!=null||r.dG_P!=null) h+=`<h3 style="font-size:13px">Binding by gorge zone</h3><div class="kv">${kv({
  "ΔG acylation site": fx(r.dG_A), "ΔG mid-gorge": fx(r.dG_M),
  "ΔG peripheral site": fx(r.dG_P), "ΔΔG P−A": fx(r.ddG)})}</div>`;
 if(r.depth){
  h+=`<h3 style="font-size:13px">Where the poses sit in this species' gorge</h3>`;
  [["chlorpyrifos_oxon","#c0392b","chlorpyrifos-oxon"],
   ["acetylcholine","#1a7f4b","acetylcholine"],
   ["TCP","#2f4b7c","TCP"]].forEach(([k,c,lab])=>{
    if(r.depth[k]) h+=`<p class="note" style="margin:6px 0 0">${lab}</p>${spark(r.depth[k],c)}`;});}
 if(r.pen) h+=`<h3 style="font-size:13px">Gorge penetration</h3><div class="kv">${kv(r.pen)}</div>`;
 d.innerHTML=h;
}
show("t0");
</script></html>"""


if __name__ == "__main__":
    main()
