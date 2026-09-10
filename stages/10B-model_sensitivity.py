#!/usr/bin/env python3
"""
Do any of the predictors explain interspecies chlorpyrifos sensitivity?

Version 2. The v1 output on this dataset showed eight structural predictors at
|r| ~ 0.6 with p_holm < 0.006, every one of which hit the Pagel lambda = 0
boundary under PGLS. v1 flagged the boundary correctly but could not say whether
that meant a real effect with phylogenetically structured residuals or a pure
taxonomic confound. It cannot: both produce lambda_residual = 0. The test that
separates them is not phylogenetic at all, it is a within-clade refit, and it is
now the centre of this script.

WHAT CHANGED, AND WHY EACH CHANGE EXISTS
----------------------------------------
0  PROVENANCE GUARD, and it is a hard stop. In the v1 run the univariate table
   carried pen_fracA_* predictors with ICCs of 0.913 / 0.893 / 0.758 while the
   penetration ICC table on disk carried pen_modefracA_* with 0.945 / 0.925 /
   0.852. Different names, different values, different run: the headline table
   predated the anchor and zone corrections in 9B-analyse_biphasic.py. Nothing in
   the pipeline noticed. Now, if a predictor table and its ICC table share no
   descriptor names, the script refuses to run and prints both name sets.

1  CLADE CONFOUND DECOMPOSITION. For every predictor: eta^2 on taxonomic class,
   and the correlation refitted inside each clade large enough to carry one. A
   predictor whose eta^2 on class EXCEEDS the response's eta^2 on class is more
   taxonomic than the thing it claims to explain, which is the arithmetic
   signature of a confound. On this dataset the structural descriptors sit at
   eta^2 0.79 to 0.96 against 0.33 for log10 LC50, and every correlation dies
   within vertebrates. That comparison is the result, and unlike PGLS it does
   not depend on a tree.

2  ICC GATING. Disattenuation by 1/sqrt(ICC) turns bottleneck_radius (ICC 0.403,
   r = -0.466) into r = -0.735, the largest number in the v1 table, produced
   entirely by dividing by the square root of a small number. Disattenuated
   values are now suppressed below --min-icc, and every predictor carries an
   explicit icc_status. Predictors with NO ICC at all, which in v1 was the whole
   geometry family including the only survivor of PGLS, are marked
   NO_NOISE_FLOOR and excluded from the reportable set unless you override.

3  SIGN CHECK. Under the binding-ranking hypothesis a more negative dG means a
   more hazardous chemical, so dG and log10 LC50 should correlate POSITIVELY.
   Five of six ligands gave a negative r in v1 and nothing said so. Affinity
   predictors now carry sign_vs_hypothesis.

4  EFFECTIVE NUMBER OF TESTS. Holm assumes independence. The structural
   descriptors have PC1 at 65 percent and Li-Ji Meff of about 7 of 11, so Holm
   is conservative and a null looks weaker than it is. Both the nominal and the
   effective-k adjustment are reported.

5  LEAVE-ONE-OUT. At n = 30 one species is 3 percent of the data. Every
   correlation is refitted n times and the range is reported.

6  ENDPOINT NOISE FLOOR. Every predictor had an ICC and the response had none.
   Supply --endpoint-variance and the reliability of log10 LC50 is estimated and
   folded into the power calculation. Without it the script says loudly that the
   response is being treated as measured without error, which it is not.

7  TREE QC. The gene tree here has 277 tips of which 27 are used, 38 percent of
   nodes below 70 percent support, 24 percent below 50, and at least one species
   present twice. Support values are now parsed and summarised, duplicate tips
   are detected and resolved by exact sequence id, and the lambda estimate
   carries a profile-likelihood interval instead of a bare point value.

8  Power docstring in v1 claimed |r| ~ 0.40 at n = 25. The function gives 0.535
   at n = 25 and 0.493 at n = 30. The printed value was right, the prose was
   wrong, and the prose is what gets read.

  module purge
  module load Python/3.11.5-GCCcore-13.2.0 SciPy-bundle/2025.07-gfbf-2025b
  python3 10B-model_sensitivity.py --endpoints ../REPORT/analysis_set.tsv \\
      --descriptors STRUCT/descriptors_species.tsv \\
      --descriptor-icc STRUCT/variance_partition.tsv \\
      --penetration BIPHASIC/penetration_descriptors.tsv \\
      --penetration-icc BIPHASIC/penetration_icc.tsv \\
      --affinity DOCKANALYSIS/affinity_species.tsv \\
      --affinity-icc DOCKANALYSIS/affinity_variance.tsv \\
      --geometry DOCKANALYSIS/nac_geometry.tsv \\
      --tree ../PHYLO/ache.raxml.support \\
      --supergroup "Vertebrata=Actinopterygii,Amphibia" \\
      --supergroup "Arthropoda=Insecta,Malacostraca,Branchiopoda" \\
      --outdir MODEL
"""

import argparse, csv, hashlib, math, os, re, sys, time
from collections import defaultdict

try:
    import numpy as np
    from scipy import stats
except ImportError:
    sys.exit("numpy+scipy required. module load SciPy-bundle/2025.07-gfbf-2025b")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_PLT = True
except ImportError:
    HAVE_PLT = False


# ------------------------------------------------------------------- io ----

def read_tsv(p):
    if not p or not os.path.exists(p):
        return []
    with open(p, newline="") as fh:
        return [{(k or "").strip(): (v.strip() if isinstance(v, str) else v)
                 for k, v in r.items()}
                for r in csv.DictReader(fh, delimiter="\t")]


def fnum(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def stamp(p):
    if not p or not os.path.exists(p):
        return ("", "", "")
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()[:12]
    return (p, time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(p))), h)


# --------------------------------------------------------------- newick ----

class Node:
    __slots__ = ("name", "length", "support", "children", "parent")
    def __init__(self):
        self.name = ""; self.length = 0.0; self.support = None
        self.children = []; self.parent = None


def parse_newick(text):
    """Newick with internal node labels captured as support values."""
    text = text.strip().rstrip(";"); i = 0

    def node():
        nonlocal i
        n = Node()
        if text[i] == "(":
            i += 1
            while True:
                c = node(); c.parent = n; n.children.append(c)
                if text[i] == ",":
                    i += 1; continue
                if text[i] == ")":
                    i += 1; break
                sys.exit(f"newick parse error near {text[max(0,i-30):i+30]!r}")
        j = i
        while i < len(text) and text[i] not in "(),:;":
            i += 1
        lab = text[j:i].strip().strip("'\"")
        if n.children:
            n.support = fnum(lab)          # RAxML / FastTree support
        else:
            n.name = lab
        if i < len(text) and text[i] == ":":
            i += 1; j = i
            while i < len(text) and (text[i].isdigit() or text[i] in ".eE+-"):
                i += 1
            n.length = fnum(text[j:i]) or 0.0
        return n
    return node()


def walk(root):
    out, st = [], [root]
    while st:
        nd = st.pop(); out.append(nd); st.extend(nd.children)
    return out


def leaves(root):
    return [n for n in walk(root) if not n.children]


def prune(root, keep):
    def rec(n):
        if not n.children:
            return n if n.name in keep else None
        kids = [k for k in (rec(c) for c in n.children) if k]
        if not kids:
            return None
        if len(kids) == 1:
            kids[0].length += n.length
            return kids[0]
        m = Node(); m.length = n.length; m.support = n.support; m.children = kids
        for k in kids:
            k.parent = m
        return m
    return rec(root)


def vcv(root, order):
    """Shared root-to-MRCA path length for every tip pair, in tree units."""
    depth, paths = {}, {}

    def rec(n, d, path):
        depth[id(n)] = d
        if not n.children:
            paths[n.name] = (path + [n], d)
        for c in n.children:
            rec(c, d + c.length, path + [n])
    rec(root, 0.0, [])
    n = len(order)
    V = np.zeros((n, n))
    for i, x in enumerate(order):
        px, dx = paths[x]
        sx = {id(z): k for k, z in enumerate(px)}
        for j, y in enumerate(order):
            py, _ = paths[y]
            last = 0.0
            for z in py:
                if id(z) in sx:
                    last = depth[id(z)]
            V[i, j] = last
        V[i, i] = dx
    return V


# ---------------------------------------------------------------- stats ----

def fisher_ci(r, n, alpha=0.05):
    if n < 4 or abs(r) >= 1:
        return (float("nan"), float("nan"))
    z = 0.5 * math.log((1 + r) / (1 - r)); se = 1 / math.sqrt(n - 3)
    zc = stats.norm.ppf(1 - alpha / 2)
    return tuple(math.tanh(v) for v in (z - zc * se, z + zc * se))


def min_detectable_r(n, alpha=0.05, power=0.80):
    if n < 5:
        return float("nan")
    za, zb = stats.norm.ppf(1 - alpha / 2), stats.norm.ppf(power)
    return math.tanh((za + zb) / math.sqrt(n - 3))


def expected_max_null(n, k):
    """Expected largest |r| among k independent tests under the null."""
    if k < 1 or n < 5:
        return float("nan")
    return math.tanh(stats.norm.ppf(1 - 0.5 / k) / math.sqrt(n - 3))


def holm(pvals, k_eff=None):
    """Holm-Bonferroni. k_eff allows the effective number of independent
    tests to be used in place of the nominal count."""
    k = len(pvals)
    keff = float(k if k_eff is None else max(1.0, min(k_eff, k)))
    idx = sorted(range(k), key=lambda i: pvals[i])
    out, prev = [0.0] * k, 0.0
    for rank, i in enumerate(idx):
        adj = min(1.0, (keff - rank) * pvals[i])
        prev = max(prev, adj); out[i] = prev
    return out


def meff_li_ji(X):
    """Li & Ji (2005) effective number of independent tests."""
    if X.shape[1] < 2 or X.shape[0] < 3:
        return float(X.shape[1])
    C = np.corrcoef(X.T)
    C = np.nan_to_num(C, nan=0.0)
    ev = np.clip(np.linalg.eigvalsh(C), 0.0, None)
    return float(sum((1.0 if e >= 1 else 0.0) + (e - math.floor(e)) for e in ev))


def eta2(y, groups):
    """Fraction of variance in y explained by a categorical grouping."""
    y = np.asarray(y, float)
    grand = y.mean()
    sst = float(((y - grand) ** 2).sum())
    if sst <= 0:
        return float("nan")
    ssb = 0.0
    for g in set(groups):
        v = y[np.array([x == g for x in groups])]
        if len(v):
            ssb += len(v) * (v.mean() - grand) ** 2
    return float(ssb / sst)


def loo_range(x, y):
    """Min and max Pearson r over all leave-one-out refits."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 6:
        return (float("nan"), float("nan"), float("nan"))
    rs, ps = [], []
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        if np.std(x[m]) == 0:
            continue
        r, p = stats.pearsonr(x[m], y[m])
        rs.append(r); ps.append(p)
    if not rs:
        return (float("nan"), float("nan"), float("nan"))
    return (min(rs), max(rs), max(ps))


def pgls(y, X, V, grid=201, fixed_lam=None):
    """GLS with Pagel's lambda. X includes the intercept.

    fixed_lam pins lambda instead of estimating it. When the predictor carries
    the same phylogenetic structure as the response, joint ML drives the
    RESIDUAL lambda to the 0 boundary, PGLS degenerates to OLS, and the
    'corrected' p comes back equal to or smaller than the uncorrected one. Both
    a causal relationship and a purely confounded one produce lambda = 0, so an
    estimated-lambda fit cannot tell them apart. Refitting at the RESPONSE's
    lambda is the cheap guard; the within-clade block is the real one.
    """
    n = len(y)
    diag = np.diag(np.diag(V)); off = V - diag

    def fit(lam):
        C = diag + lam * off + 1e-9 * np.eye(n)
        try:
            L = np.linalg.cholesky(C)
        except np.linalg.LinAlgError:
            return None
        Ci = np.linalg.inv(C)
        XtCi = X.T @ Ci
        try:
            cov = np.linalg.inv(XtCi @ X)
        except np.linalg.LinAlgError:
            return None
        b = cov @ (XtCi @ y)
        r = y - X @ b
        s2 = float(r.T @ Ci @ r) / (n - X.shape[1])
        s2ml = float(r.T @ Ci @ r) / n
        ll = -0.5 * (n * math.log(2 * math.pi * s2ml)
                     + 2 * float(np.sum(np.log(np.diag(L)))) + n)
        return b, s2 * cov, ll

    profile, best = [], None
    grid_vals = ([float(fixed_lam)] if fixed_lam is not None
                 else np.linspace(0, 1, grid))
    for lam in grid_vals:
        f = fit(float(lam))
        if f is None:
            continue
        profile.append((float(lam), f[2]))
        if best is None or f[2] > best[1][2]:
            best = (float(lam), f)
    if best is None:
        return None
    lam, (b, cov, ll) = best
    se = math.sqrt(max(cov[-1, -1], 0.0))
    t = b[-1] / se if se > 0 else float("nan")
    p = 2 * stats.t.sf(abs(t), n - X.shape[1]) if se > 0 else float("nan")
    lo = hi = lam
    if len(profile) > 1:
        ok = [l for l, v in profile if v >= ll - 1.92]
        if ok:
            lo, hi = min(ok), max(ok)
    return {"lambda": lam, "lambda_lo": lo, "lambda_hi": hi,
            "beta": float(b[-1]), "se": se, "t": float(t), "p": float(p),
            "logLik": ll, "profile": profile}


# ------------------------------------------------------------------ main ---

FAMILY_ICC_ARG = {"structural": "--descriptor-icc",
                  "penetration": "--penetration-icc",
                  "affinity": "--affinity-icc",
                  "geometry": "(none supplied)"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--endpoints", required=True)
    ap.add_argument("--descriptors")
    ap.add_argument("--descriptor-icc")
    ap.add_argument("--penetration")
    ap.add_argument("--penetration-icc")
    ap.add_argument("--affinity")
    ap.add_argument("--affinity-icc")
    ap.add_argument("--geometry")
    ap.add_argument("--geometry-icc",
                    help="per-model / per-seed ICC for the near-attack "
                         "geometry columns. Without it the whole geometry "
                         "family has no noise floor and is excluded from the "
                         "reportable set.")
    ap.add_argument("--endpoint-variance",
                    help="TSV with columns species (tag or ECOTOX name), "
                         "n_records, sd_log10. Gives the response a "
                         "reliability estimate and folds it into power.")
    ap.add_argument("--tree")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-icc", type=float, default=0.75,
                    help="ICC below which disattenuated r is suppressed and "
                         "the predictor is not reportable")
    ap.add_argument("--min-clade-n", type=int, default=8,
                    help="smallest clade in which a within-clade correlation "
                         "is attempted")
    ap.add_argument("--supergroup", action="append", default=[],
                    help='"Name=ClassA,ClassB". Repeatable.')
    ap.add_argument("--allow-no-icc", action="store_true",
                    help="let predictors with no noise floor into the "
                         "reportable set. Off by default.")
    ap.add_argument("--no-provenance-stop", action="store_true",
                    help="downgrade the stale-input check to a warning")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    # ---- 0. provenance -----------------------------------------------------
    print("== 0. provenance ==")
    inputs = [("endpoints", a.endpoints), ("descriptors", a.descriptors),
              ("descriptor_icc", a.descriptor_icc),
              ("penetration", a.penetration),
              ("penetration_icc", a.penetration_icc),
              ("affinity", a.affinity), ("affinity_icc", a.affinity_icc),
              ("geometry", a.geometry), ("geometry_icc", a.geometry_icc),
              ("tree", a.tree)]
    prov = []
    for label, p in inputs:
        path, mt, h = stamp(p)
        if path:
            prov.append({"role": label, "path": path, "modified": mt,
                         "sha256_12": h})
            print(f"  {label:16s} {mt}  {h}  {path}")
        elif p:
            print(f"  {label:16s} MISSING  {p}")
    with open(os.path.join(a.outdir, "provenance.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=["role", "path", "modified", "sha256_12"])
        w.writeheader(); w.writerows(prov)

    # ---- endpoints ---------------------------------------------------------
    eps = {}
    for r in read_tsv(a.endpoints):
        sid = (r.get("Sequence_id") or r.get("sequence_id") or "").strip()
        y = fnum(r.get("log10_LC50"))
        if not sid or y is None:
            continue
        eps[sid] = {"sequence_id": sid, "tag": sid.split("__")[0],
                    "species": (r.get("Species_representative")
                                or r.get("species") or sid.split("__")[0]),
                    "cls": r.get("Class", ""), "y": y,
                    "n_collapsed": int(fnum(r.get("n_species_collapsed")) or 1),
                    "members": r.get("member_species", "")}
    if not eps:
        sys.exit("no endpoints parsed; check Sequence_id / log10_LC50")
    by_tag = {v["tag"]: v for v in eps.values()}

    # supergroups
    sg = {}
    for spec in a.supergroup:
        if "=" not in spec:
            sys.exit(f'--supergroup needs Name=ClassA,ClassB, got {spec!r}')
        nm, cls = spec.split("=", 1)
        for c in cls.split(","):
            sg.setdefault(c.strip(), []).append(nm.strip())

    # ---- predictors --------------------------------------------------------
    icc, icc_src = {}, {}

    def load_icc(path, keycol, prefix="", family=""):
        got = {}
        for r in read_tsv(path):
            k = (r.get(keycol) or "").strip()
            v = fnum(r.get("ICC"))
            if k and v is not None:
                got[prefix + k] = v
        for k, v in got.items():
            icc[k] = v; icc_src[k] = family
        return set(got)

    icc_names = {}
    icc_names["structural"] = load_icc(a.descriptor_icc, "descriptor",
                                       family="structural")
    icc_names["penetration"] = load_icc(a.penetration_icc, "descriptor",
                                        family="penetration")
    icc_names["affinity"] = load_icc(a.affinity_icc, "ligand", "aff_",
                                     family="affinity")
    icc_names["geometry"] = load_icc(a.geometry_icc, "descriptor", "geom_",
                                     family="geometry")

    pred = defaultdict(dict)
    fam = {}

    def take_table(path, family, skip=(), prefix=""):
        names = set()
        for r in read_tsv(path):
            tag = (r.get("species") or "").strip()
            if not tag:
                continue
            for k, v in r.items():
                if k in skip or k.endswith("_sd") or k in ("species",):
                    continue
                fv = fnum(v)
                if fv is None:
                    continue
                nm = prefix + k
                pred[tag][nm] = fv; fam[nm] = family; names.add(nm)
        return names

    tab_names = {}
    tab_names["structural"] = take_table(a.descriptors, "structural",
                                         skip=("sequence_id", "n_models"))
    tab_names["penetration"] = take_table(a.penetration, "penetration")

    aff_names = set()
    for r in read_tsv(a.affinity):
        tag = (r.get("species") or "").strip()
        lig = (r.get("ligand") or "").strip()
        fv = fnum(r.get("affinity_mean"))
        if tag and lig and fv is not None:
            nm = "aff_" + lig
            pred[tag][nm] = fv; fam[nm] = "affinity"; aff_names.add(nm)
    tab_names["affinity"] = aff_names

    geom_names = set()
    GEOM_COLS = ("d_SerOG_P", "inline_deviation", "oxyanion_min_dist",
                 "n_oxyanion_hbond", "n_contact_residues")
    for r in read_tsv(a.geometry):
        tag = (r.get("species") or "").strip()
        for k in GEOM_COLS:
            fv = fnum(r.get(k))
            if tag and fv is not None:
                nm = "geom_" + k
                pred[tag][nm] = fv; fam[nm] = "geometry"; geom_names.add(nm)
    tab_names["geometry"] = geom_names

    # ---- 0b. stale-input guard --------------------------------------------
    stale = []
    for family, names in tab_names.items():
        iccs = icc_names.get(family, set())
        if not names or not iccs:
            continue
        inter = {n for n in names if n in iccs}
        if not inter:
            stale.append((family, sorted(names)[:6], sorted(iccs)[:6]))
    if stale:
        print("\n*** STALE INPUT DETECTED ***")
        for family, n1, n2 in stale:
            print(f"  family '{family}': the predictor table and the ICC table "
                  f"({FAMILY_ICC_ARG.get(family,'?')}) share NO descriptor "
                  f"names.")
            print(f"    predictor table : {n1}")
            print(f"    ICC table       : {n2}")
        print("  These two files came from different runs. Any ICC-gated or "
              "disattenuated result would be attached to the wrong predictor.")
        if not a.no_provenance_stop:
            sys.exit("  Regenerate both from the same run, or pass "
                     "--no-provenance-stop if you understand the consequence.")
        print("  CONTINUING because --no-provenance-stop was given.\n")

    # ---- join --------------------------------------------------------------
    matched = sorted(t for t in by_tag if t in pred)
    missing = sorted(t for t in by_tag if t not in pred)
    print(f"\n== 1. join ==")
    print(f"  endpoint species {len(by_tag)}; with predictors {len(matched)}")
    if missing:
        print("  endpoint species with NO predictors (excluded):")
        for t in missing:
            print(f"    {by_tag[t]['species']:32s} tag={t}")
    if len(matched) < 6:
        sys.exit("too few joined species to model")

    with open(os.path.join(a.outdir, "join_report.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sequence_id", "species", "tag", "class", "log10_LC50",
                    "joined", "n_predictors", "n_species_collapsed"])
        for t, e in sorted(by_tag.items()):
            w.writerow([e["sequence_id"], e["species"], t, e["cls"], e["y"],
                        t in pred, len(pred.get(t, {})), e["n_collapsed"]])

    names = sorted({k for t in matched for k in pred[t]})
    n_full = len(matched)
    yv_all = np.array([by_tag[t]["y"] for t in matched])
    cls_all = [by_tag[t]["cls"] for t in matched]

    # ---- 2. the response's own noise floor --------------------------------
    print(f"\n== 2. the response ==")
    print(f"  log10 LC50 {yv_all.min():.2f} to {yv_all.max():.2f} "
          f"(span {yv_all.max()-yv_all.min():.2f} log units), SD "
          f"{yv_all.std(ddof=1):.2f}")
    from collections import Counter
    for c, k in Counter(cls_all).most_common():
        print(f"    {c}: {k}")
    coll = [t for t in matched if by_tag[t]["n_collapsed"] > 1]
    if coll:
        print(f"  {len(coll)} endpoints are means over more than one ECOTOX "
              f"species:")
        for t in coll:
            print(f"    {by_tag[t]['species']:26s} n={by_tag[t]['n_collapsed']}"
                  f"  {by_tag[t]['members']}")

    rel_y = None
    if a.endpoint_variance:
        ev = {}
        for r in read_tsv(a.endpoint_variance):
            k = (r.get("species") or "").strip().replace(" ", "_")
            nrec = fnum(r.get("n_records")) or 1.0
            sd = fnum(r.get("sd_log10"))
            if k and sd is not None:
                ev[k] = (nrec, sd)
        use = [(by_tag[t]["y"], *ev[t]) for t in matched if t in ev]
        if len(use) >= 6:
            ys = np.array([u[0] for u in use])
            sem2 = np.array([(u[2] ** 2) / max(u[1], 1.0) for u in use])
            vb = max(0.0, float(ys.var(ddof=1)) - float(sem2.mean()))
            rel_y = vb / (vb + float(sem2.mean())) if (vb + sem2.mean()) > 0 else 0.0
            print(f"  endpoint reliability from {len(use)} species with "
                  f"replicate records: R_y = {rel_y:.3f}")
            print(f"    mean within-species SD {math.sqrt(float(sem2.mean())*1):.2f} "
                  f"log units on the mean")
        else:
            print("  --endpoint-variance supplied but fewer than 6 species "
                  "matched; response reliability not estimated")
    if rel_y is None:
        print("  *** the response has NO noise floor in this run ***")
        print("  log10 LC50 is a geometric mean over ECOTOX records that differ "
              "in duration, exposure system, nominal versus measured "
              "concentration and life stage. Treating it as measured without "
              "error overstates the power figures below. Supply "
              "--endpoint-variance.")

    # ---- 3. power ----------------------------------------------------------
    mdr = min_detectable_r(n_full)
    print(f"\n== 3. power ==")
    print(f"  n = {n_full}; minimum |r| detectable at alpha 0.05, power 0.80: "
          f"{mdr:.3f}")
    fam_names = defaultdict(list)
    for nm in names:
        fam_names[fam.get(nm, "?")].append(nm)
    keff = {}
    for f, ns in sorted(fam_names.items()):
        M = []
        for t in matched:
            row = [pred[t].get(n) for n in ns]
            if all(v is not None for v in row):
                M.append(row)
        X = np.array(M, float)
        k = meff_li_ji(X) if X.shape[0] >= 3 else float(len(ns))
        keff[f] = k
        print(f"    family {f:12s} {len(ns):2d} predictors, Li-Ji effective "
              f"{k:.1f}   expected max |r| under the null "
              f"{expected_max_null(n_full, max(1, round(k))):.2f}")
    if rel_y is not None:
        print(f"  with R_y = {rel_y:.2f}, an observed |r| of {mdr:.2f} "
              f"corresponds to a TRUE |r| of about "
              f"{min(1.0, mdr/math.sqrt(rel_y)):.2f} once endpoint error is "
              f"accounted for.")

    # ---- 4. tree -----------------------------------------------------------
    tree_res, V, order = None, None, None
    if a.tree and os.path.exists(a.tree):
        root = parse_newick(open(a.tree).read())
        allt = leaves(root)
        sups = [n.support for n in walk(root)
                if n.children and n.support is not None]
        print(f"\n== 4. tree ==")
        print(f"  {len(allt)} tips in file; {len(by_tag)} analysis-set species")
        if sups:
            s = np.array(sups, float)
            if s.max() <= 1.0:
                s = s * 100.0
            print(f"  node support: median {np.median(s):.0f}, "
                  f"{100*np.mean(s < 70):.0f}% below 70, "
                  f"{100*np.mean(s < 50):.0f}% below 50")
            if np.mean(s < 70) > 0.25:
                print("  A quarter or more of the topology is weakly "
                      "supported. Branch lengths from this tree carry that "
                      "uncertainty into the VCV, and lambda inherits it. "
                      "Report lambda with its profile interval, not as a "
                      "point estimate.")

        sid_of = {t: by_tag[t]["sequence_id"] for t in matched}
        exact = {v: k for k, v in sid_of.items()}
        # exact sequence ids first; tags only where unambiguous
        tip2tag, dupes = {}, defaultdict(list)
        for lf in allt:
            if lf.name in exact:
                tip2tag[lf.name] = exact[lf.name]
            base = lf.name.split("__")[0]
            if base in by_tag:
                dupes[base].append(lf.name)
        for tg, tips in dupes.items():
            if tg in tip2tag.values():
                continue
            if len(tips) == 1:
                tip2tag[tips[0]] = tg
            else:
                print(f"  DUPLICATE TIPS for {tg}: {tips}")
                print(f"    none is the analysis-set sequence id "
                      f"({sid_of.get(tg,'?')}); skipped rather than guessed.")
        multi = {tg: v for tg, v in dupes.items() if len(v) > 1}
        if multi:
            per = sorted(len(v) for v in multi.values())
            print(f"  {len(multi)} of {len(by_tag)} analysis-set species have "
                  f"more than one tip on this tree (median "
                  f"{per[len(per)//2]}, max {per[-1]} sequences per species). "
                  f"This is a gene tree of a multi-gene family, not a species "
                  f"phylogeny. One tip per species is selected by exact "
                  f"analysis-set sequence id; the choice of tip changes the "
                  f"VCV, so it is an assumption, not a lookup.")

        want = set(tip2tag)
        looks_gene = any("__" in t for t in want)
        pr = prune(root, want)
        if pr is not None:
            for nd in leaves(pr):
                nd.name = tip2tag.get(nd.name, nd.name)
            order = [nd.name for nd in leaves(pr)]
            order = [t for t in order if t in by_tag]
            print(f"  on tree after pruning: {len(order)} of {len(matched)}")
            dropped = [t for t in matched if t not in set(order)]
            if dropped:
                print(f"  NOT on the tree ({len(dropped)}), excluded from PGLS "
                      f"only: {', '.join(dropped)}")
            if len(order) >= 6:
                V = vcv(pr, order)
                yv = np.array([by_tag[s]["y"] for s in order])
                r0 = pgls(yv, np.ones((len(order), 1)), V)
                if r0:
                    Cn = np.diag(np.diag(V)); Cin = np.linalg.inv(Cn)
                    Xi = np.ones((len(order), 1))
                    b = np.linalg.inv(Xi.T @ Cin @ Xi) @ (Xi.T @ Cin @ yv)
                    res = yv - Xi @ b
                    s2 = float(res.T @ Cin @ res) / len(yv)
                    sign, logdet = np.linalg.slogdet(Cn)
                    ll0 = -0.5 * (len(yv) * math.log(2 * math.pi * s2)
                                  + float(logdet) + len(yv))
                    tree_res = {"lambda": r0["lambda"],
                                "lambda_lo": r0["lambda_lo"],
                                "lambda_hi": r0["lambda_hi"],
                                "logLik": r0["logLik"], "logLik_lambda0": ll0,
                                "LR": 2 * (r0["logLik"] - ll0),
                                "n": len(order),
                                "tree_type": "gene" if looks_gene else "species"}
                    print(f"\n  phylogenetic signal in log10 LC50: "
                          f"lambda = {tree_res['lambda']:.3f} "
                          f"[{tree_res['lambda_lo']:.2f}, "
                          f"{tree_res['lambda_hi']:.2f}]  "
                          f"LR vs 0 = {tree_res['LR']:.2f} (1 df)")
                    if looks_gene:
                        print("  CAVEAT: tips are sequence ids, so this is the "
                              "GENE tree. Branch lengths are substitutions per "
                              "site, not time, and every structural predictor "
                              "is a function of the same sequences. "
                              "Diagnostic only; rerun with the reconciled "
                              "TimeTree.")
        if V is None:
            print("\n  *** TREE SUPPLIED BUT PGLS DID NOT RUN ***")
            print("  No phylogenetic correction was applied. Every PGLS column "
                  "below is empty.")

    # ---- 5. univariate + clade decomposition ------------------------------
    eta_y = eta2(yv_all, cls_all)
    clade_defs = []
    counts = Counter(cls_all)
    for c, k in counts.items():
        if k >= a.min_clade_n:
            clade_defs.append((c, [t for t in matched if by_tag[t]["cls"] == c]))
    sup_members = defaultdict(list)
    for t in matched:
        for nm in sg.get(by_tag[t]["cls"], []):
            sup_members[nm].append(t)
    for nm, ts in sup_members.items():
        if len(ts) >= a.min_clade_n:
            clade_defs.append((nm, ts))
    clade_defs.sort(key=lambda kv: -len(kv[1]))

    rows = []
    for nm in names:
        use_t = [t for t in matched if nm in pred[t]]
        if len(use_t) < 6:
            continue
        x = np.array([pred[t][nm] for t in use_t], float)
        y = np.array([by_tag[t]["y"] for t in use_t], float)
        if np.std(x) == 0:
            continue
        pr_, pp = stats.pearsonr(x, y)
        sr, _ = stats.spearmanr(x, y)
        lo, hi = fisher_ci(pr_, len(x))
        ic = icc.get(nm)
        f = fam.get(nm, "?")

        if ic is None:
            icc_status = "NO_NOISE_FLOOR"
        elif ic >= a.min_icc:
            icc_status = "ok"
        elif ic >= 0.5:
            icc_status = "marginal"
        else:
            icc_status = "mostly_noise"

        row = {"predictor": nm, "family": f, "n": len(x),
               "pearson_r": round(float(pr_), 3), "p_raw": float(pp),
               "ci_lo": round(lo, 3), "ci_hi": round(hi, 3),
               "spearman_rho": round(float(sr), 3),
               "ICC": round(ic, 3) if ic is not None else "",
               "icc_status": icc_status,
               "detectable": "yes" if abs(pr_) >= mdr else "below_power"}

        # disattenuation only where the noise floor supports it
        if ic is not None and ic >= a.min_icc:
            row["r_disattenuated"] = round(float(pr_) / math.sqrt(ic), 3)
        else:
            row["r_disattenuated"] = ""

        # sign against the binding-ranking hypothesis
        if f == "affinity":
            row["sign_vs_hypothesis"] = ("as_expected" if pr_ > 0
                                         else "REVERSED")
        else:
            row["sign_vs_hypothesis"] = ""

        # leave-one-out
        lmin, lmax, lpmax = loo_range(x, y)
        row.update(loo_r_min=round(lmin, 3), loo_r_max=round(lmax, 3),
                   loo_p_max=round(lpmax, 4),
                   loo_sign_stable=bool(lmin * lmax > 0))

        # clade confound
        cls_u = [by_tag[t]["cls"] for t in use_t]
        e_x = eta2(x, cls_u)
        row["eta2_class"] = round(e_x, 3)
        row["eta2_class_ratio"] = (round(e_x / eta_y, 2)
                                   if eta_y and eta_y > 0 else "")
        row["more_taxonomic_than_endpoint"] = bool(eta_y and e_x > eta_y)
        # bounded or censored predictors: a fraction stacked at the observed
        # extreme makes a Pearson r on them close to meaningless, and the
        # within-group variance used by ICC collapses, inflating the ICC.
        frac_ext = float(max((x == x.max()).mean(), (x == x.min()).mean()))
        row["frac_at_extreme"] = round(frac_ext, 3)
        row["censored"] = bool(frac_ext >= 0.25)

        best_within, n_within, best_p = 0.0, 0, 1.0
        for cname, members in clade_defs:
            sub = [t for t in members if nm in pred[t]]
            if len(sub) < 6:
                continue
            xs = np.array([pred[t][nm] for t in sub], float)
            ys = np.array([by_tag[t]["y"] for t in sub], float)
            if np.std(xs) == 0:
                continue
            rr, ppv = stats.pearsonr(xs, ys)
            row[f"r_in_{cname}"] = round(float(rr), 3)
            row[f"p_in_{cname}"] = round(float(ppv), 4)
            row[f"n_in_{cname}"] = len(sub)
            if abs(rr) > abs(best_within):
                best_within, n_within, best_p = float(rr), len(sub), float(ppv)
        n_cl = sum(1 for cname, _ in clade_defs if f"r_in_{cname}" in row)
        row["n_clades_tested"] = n_cl
        row["max_abs_r_within_clade"] = round(abs(best_within), 3)
        row["n_in_clade_of_max_r"] = n_within
        row["p_in_clade_of_max_r"] = round(best_p, 4)
        mdr_w = min_detectable_r(n_within) if n_within >= 6 else float("nan")
        row["mdr_within_clade"] = (round(mdr_w, 3)
                                   if math.isfinite(mdr_w) else "")
        # the within-clade r is a maximum over n_clades_tested groups, so the
        # threshold has to be corrected for that or it is a best-of-k result
        # reported as a single test.
        surv = (n_cl > 0
                and abs(best_within) >= (mdr_w if math.isfinite(mdr_w) else 1.0)
                and best_p < 0.05 / max(1, n_cl))
        if not clade_defs:
            row["clade_verdict"] = "not_tested"
        elif surv and row["censored"]:
            row["clade_verdict"] = "within_clade_but_CENSORED"
        elif surv:
            row["clade_verdict"] = "survives_within_clade"
        elif row["more_taxonomic_than_endpoint"]:
            row["clade_verdict"] = "CLASS_CONFOUNDED"
        else:
            row["clade_verdict"] = "underpowered_within_clade"

        # PGLS
        if V is not None:
            useo = [s for s in order if nm in pred[s]]
            if len(useo) >= 6:
                idx = [order.index(s) for s in useo]
                Vs = V[np.ix_(idx, idx)]
                yo = np.array([by_tag[s]["y"] for s in useo])
                xo = np.array([pred[s][nm] for s in useo])
                Xd = np.column_stack([np.ones(len(useo)), xo])
                g = pgls(yo, Xd, Vs)
                if g:
                    row.update(pgls_n=len(useo),
                               pgls_lambda=round(g["lambda"], 3),
                               pgls_beta=round(g["beta"], 4),
                               pgls_se=round(g["se"], 4),
                               pgls_p=g["p"],
                               pgls_lambda_at_boundary=bool(
                                   g["lambda"] <= 0.001 or g["lambda"] >= 0.999))
                    if tree_res is not None:
                        gf = pgls(yo, Xd, Vs, fixed_lam=tree_res["lambda"])
                        if gf:
                            row.update(pgls_fixlam=round(tree_res["lambda"], 3),
                                       pgls_fix_beta=round(gf["beta"], 4),
                                       pgls_fix_p=gf["p"])
        rows.append(row)

    if not rows:
        sys.exit("no predictor had at least 6 joined species")

    # ---- multiple testing, within family, at nominal and effective k ------
    for f in sorted({r["family"] for r in rows}):
        grp = [r for r in rows if r["family"] == f]
        k = keff.get(f)
        for key, out_nom, out_eff in (
                ("p_raw", "p_holm", "p_holm_eff"),
                ("pgls_p", "pgls_p_holm", None),
                ("pgls_fix_p", "pgls_fix_p_holm", None)):
            sub = [r for r in grp if isinstance(r.get(key), float)]
            if not sub:
                continue
            ps = [r[key] for r in sub]
            for r, q in zip(sub, holm(ps)):
                r[out_nom] = q
            if out_eff:
                for r, q in zip(sub, holm(ps, k_eff=k)):
                    r[out_eff] = q
    for r in rows:
        for key in ("p_raw", "p_holm", "p_holm_eff", "pgls_p", "pgls_p_holm",
                    "pgls_fix_p", "pgls_fix_p_holm"):
            if isinstance(r.get(key), float):
                r[key] = round(r[key], 4)

    rows.sort(key=lambda r: (r["family"], -abs(r["pearson_r"])))
    clade_cols = []
    for cname, _ in clade_defs:
        clade_cols += [f"n_in_{cname}", f"r_in_{cname}", f"p_in_{cname}"]
    cols = ["predictor", "family", "n", "pearson_r", "ci_lo", "ci_hi",
            "spearman_rho", "p_raw", "p_holm", "p_holm_eff",
            "ICC", "icc_status", "r_disattenuated", "sign_vs_hypothesis",
            "detectable", "loo_r_min", "loo_r_max", "loo_p_max",
            "loo_sign_stable",
            "eta2_class", "eta2_class_ratio", "more_taxonomic_than_endpoint",
            "frac_at_extreme", "censored", "n_clades_tested",
            "max_abs_r_within_clade", "n_in_clade_of_max_r",
            "p_in_clade_of_max_r", "mdr_within_clade",
            "clade_verdict"] + clade_cols + [
            "pgls_n", "pgls_lambda", "pgls_lambda_at_boundary", "pgls_beta",
            "pgls_se", "pgls_p", "pgls_p_holm", "pgls_fixlam",
            "pgls_fix_beta", "pgls_fix_p", "pgls_fix_p_holm"]
    with open(os.path.join(a.outdir, "univariate_results.tsv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                           fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    # ---- 6. report ---------------------------------------------------------
    print(f"\n== 5. predictors vs log10 LC50 ==")
    print(f"{'predictor':28s} {'fam':11s} {'r':>7s} {'p_holm':>7s} "
          f"{'ICC':>6s} {'eta2cl':>7s} {'r_in_clade':>11s} {'verdict':>22s}")
    for r in rows:
        icv = r["ICC"] if r["ICC"] != "" else float("nan")
        print(f"{r['predictor']:28s} {r['family']:11s} "
              f"{r['pearson_r']:+7.3f} {r.get('p_holm',1.0):7.3f} "
              f"{icv if isinstance(icv,float) else 0:6.3f} "
              f"{r['eta2_class']:7.2f} "
              f"{r['max_abs_r_within_clade']:11.3f} "
              f"{r['clade_verdict']:>22s}")

    print(f"\n== 6. clade confound ==")
    print(f"  eta^2 of log10 LC50 on class: {eta_y:.3f}")
    for cname, members in clade_defs:
        ys = np.array([by_tag[t]["y"] for t in members])
        print(f"    clade {cname:16s} n={len(members):2d}  "
              f"log10 LC50 {ys.min():+.2f} to {ys.max():+.2f}  "
              f"detectable |r| {min_detectable_r(len(members)):.2f}")
    worse = [r for r in rows if r["more_taxonomic_than_endpoint"]]
    print(f"\n  {len(worse)} of {len(rows)} predictors are MORE taxonomic than "
          f"the endpoint they claim to explain (eta^2 on class exceeds "
          f"{eta_y:.2f}):")
    for r in sorted(worse, key=lambda x: -x["eta2_class"])[:14]:
        print(f"    {r['predictor']:28s} eta^2 {r['eta2_class']:.2f} "
              f"({r['eta2_class_ratio']}x the endpoint)  "
              f"r all {r['pearson_r']:+.2f} -> within clade "
              f"{r['max_abs_r_within_clade']:.2f}")
    print("\n  A predictor that partitions on class more sharply than the "
          "response does is separating taxa, not explaining sensitivity. This "
          "comparison needs no tree and no lambda.")

    surv = [r for r in rows if r["clade_verdict"] == "survives_within_clade"]
    print(f"\n  predictors surviving inside a clade: {len(surv)}")
    for r in surv:
        print(f"    {r['predictor']} r={r['max_abs_r_within_clade']:.3f} "
              f"in n={r['n_in_clade_of_max_r']} (detectable "
              f"{r['mdr_within_clade']}, best of {r['n_clades_tested']} "
              f"clades)")

    # ---- 7. verdict --------------------------------------------------------
    print(f"\n== 7. what is reportable ==")
    def reportable(r):
        if r["icc_status"] == "NO_NOISE_FLOOR" and not a.allow_no_icc:
            return False
        if r["icc_status"] == "mostly_noise":
            return False
        if r["censored"]:
            return False
        if not r["loo_sign_stable"]:
            return False
        return r.get("p_holm", 1.0) < 0.05
    rep = [r for r in rows if reportable(r)]
    noicc = [r for r in rows if r["icc_status"] == "NO_NOISE_FLOOR"]
    if noicc:
        print(f"  {len(noicc)} predictors have NO noise floor "
              f"({', '.join(sorted({r['family'] for r in noicc}))} family): "
              f"{', '.join(r['predictor'] for r in noicc)}")
        print("  They are computed from a single pose or a single model with "
              "no replicate. Give them the same replicate treatment as the "
              "structural descriptors or do not report them.")
    cen = [r for r in rows if r["censored"]]
    if cen:
        print(f"  {len(cen)} predictors are censored (a quarter or more of "
              f"species sit at the observed extreme): "
              f"{', '.join(r['predictor'] for r in cen)}")
        print("  A Pearson r on a stacked variable is not interpretable, and "
              "the within-species variance that the ICC divides by collapses "
              "at the bound, so their ICCs are inflated too. Fix the "
              "definition rather than the statistics.")
    print(f"  passing Holm, with an adequate noise floor and a stable "
          f"leave-one-out sign: {len(rep)}")
    for r in rep:
        print(f"    {r['predictor']:28s} r={r['pearson_r']:+.3f} "
              f"p_holm={r.get('p_holm'):.4f} ICC={r['ICC']} "
              f"clade_verdict={r['clade_verdict']}")
    rep_clade = [r for r in rep if r["clade_verdict"] == "survives_within_clade"]
    nb = [r for r in rows if r.get("pgls_lambda_at_boundary")]
    if nb:
        print(f"\n  {len(nb)} predictors hit the lambda boundary under PGLS. "
              f"For these the 'corrected' p is not a correction. Read the "
              f"clade_verdict column instead: it distinguishes the two cases "
              f"that lambda = 0 cannot.")
    if not rep_clade:
        print(f"\n  NO predictor survives inside a clade. State the result as "
              f"a bound: across {n_full} species spanning "
              f"{yv_all.max()-yv_all.min():.1f} log units, the structural and "
              f"docking predictors track taxonomic class "
              f"(median eta^2 "
              f"{np.median([r['eta2_class'] for r in rows]):.2f}) rather than "
              f"sensitivity (eta^2 {eta_y:.2f}), and within the largest clade "
              f"nothing reaches the detectable threshold. That is a "
              f"quantitative negative result about target-site structure, not "
              f"a failed analysis.")
    aff_rev = [r for r in rows
               if r["sign_vs_hypothesis"] == "REVERSED"]
    if aff_rev:
        print(f"\n  {len(aff_rev)} of "
              f"{len([r for r in rows if r['family']=='affinity'])} affinity "
              f"predictors have the sign REVERSED against the binding-ranking "
              f"hypothesis (tighter predicted binding going with LOWER "
              f"sensitivity): "
              f"{', '.join(r['predictor'] for r in aff_rev)}")
        print("  Under BARMIE logic dG and log10 LC50 should correlate "
              "positively. Explain this rather than omitting it; the usual "
              "cause is that the score tracks pocket volume and pocket volume "
              "tracks taxon.")

    # ---- 8. figures --------------------------------------------------------
    if HAVE_PLT and rows:
        colr = {"structural": "#2c7fb8", "affinity": "#e8703a",
                "penetration": "#7b5aa6", "geometry": "#4a9b5c"}
        missing_c = sorted({r["family"] for r in rows} - set(colr))
        if missing_c:
            print(f"  NOTE: families with no assigned colour: {missing_c}")

        fig, ax = plt.subplots(figsize=(8.0, 0.32 * len(rows) + 2))
        for i, r in enumerate(rows):
            c = colr.get(r["family"], "0.5")
            ax.plot([r["ci_lo"], r["ci_hi"]], [i, i], color=c, lw=2, alpha=.75)
            ax.scatter([r["pearson_r"]], [i], s=30, color=c, zorder=3)
            if r["max_abs_r_within_clade"]:
                sgn = 1 if r["pearson_r"] >= 0 else -1
                ax.scatter([sgn * r["max_abs_r_within_clade"]], [i], s=26,
                           facecolors="none", edgecolors=c, zorder=3)
        ax.axvline(0, color="0.3", lw=0.8)
        for s in (-1, 1):
            ax.axvline(s * mdr, color="#c0392b", ls="--", lw=0.9)
        ax.text(mdr, -0.9, f" detectable |r|={mdr:.2f}", fontsize=7,
                color="#c0392b", va="top")
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([r["predictor"] for r in rows], fontsize=7)
        ax.invert_yaxis(); ax.set_xlabel("Pearson r vs log10 LC50")
        ax.set_title(f"Predictors of acute chlorpyrifos sensitivity "
                     f"(n={n_full}); open circles = best within-clade |r|",
                     fontsize=9)
        present = [f for f in colr if any(r["family"] == f for r in rows)]
        ax.legend(handles=[plt.Line2D([], [], color=colr[k], lw=2, label=k)
                           for k in present],
                  fontsize=7, frameon=False, loc="lower right")
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "fig_forest.png"), dpi=200)
        fig.savefig(os.path.join(a.outdir, "fig_forest.svg"))
        plt.close(fig)

        # the confound figure
        fig, ax = plt.subplots(figsize=(6.4, 5.0))
        for r in rows:
            c = colr.get(r["family"], "0.5")
            ax.scatter(r["eta2_class"], abs(r["pearson_r"]), s=34, color=c,
                       zorder=3)
        ax.axvline(eta_y, color="#c0392b", ls="--", lw=1.0)
        ax.text(eta_y, ax.get_ylim()[1], f" eta$^2$ of endpoint = {eta_y:.2f}",
                fontsize=7, color="#c0392b", va="top")
        ax.axhline(mdr, color="0.4", ls=":", lw=0.9)
        ax.set_xlabel("fraction of predictor variance explained by "
                      "taxonomic class")
        ax.set_ylabel("|r| with log10 LC50")
        ax.set_title("Predictors to the right of the red line are more\n"
                     "taxonomic than the endpoint they explain", fontsize=9)
        ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=colr[k],
                                      label=k) for k in present],
                  fontsize=7, frameon=False, loc="lower right")
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "fig_clade_confound.png"), dpi=200)
        fig.savefig(os.path.join(a.outdir, "fig_clade_confound.svg"))
        plt.close(fig)

    if tree_res:
        tr = {k: v for k, v in tree_res.items() if k != "profile"}
        with open(os.path.join(a.outdir, "phylo_signal.tsv"), "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                               fieldnames=list(tr))
            w.writeheader(); w.writerow(tr)

    print(f"\n  {a.outdir}/univariate_results.tsv")
    print(f"  {a.outdir}/join_report.tsv")
    print(f"  {a.outdir}/provenance.tsv")
    if HAVE_PLT:
        print(f"  {a.outdir}/fig_forest.png")
        print(f"  {a.outdir}/fig_clade_confound.png")


if __name__ == "__main__":
    main()
