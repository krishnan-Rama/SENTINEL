#!/usr/bin/env python3
"""
Phylogenetics, active-site analysis, SSD, and a figure at each step.

Steps, each writing a figure:
  1  prune the AChE tree to the analysis species                fig1_tree
  2  active-site residue matrix across species                  fig2_active_site
  3  model coverage per sequence (domain completeness)          fig3_coverage
  4  Pagel's lambda for log10 LC50 on the pruned tree           fig4_lambda
  5  species sensitivity distribution and HC5                   fig5_ssd
  6  LC50 by class and by clade                                 fig6_by_class

Read step 4 before you invest in docking. If lambda is near 1, almost all the
interspecies variance in acute sensitivity is explained by shared ancestry, and
whatever a structural predictor can add is confined to the residual. That is a
publishable result either way, and it costs a minute here against months of
docking.

One caveat you must carry: the tree used here is the AChE GENE tree, not a dated
species tree. A gene tree conflates sequence divergence with time and can be
distorted by paralogy and rate variation, so the lambda estimate is a first-pass
diagnostic, not the final number. Pass --tree with a dated species tree
(TimeTree or similar) for the version that goes in the paper.

  module purge
  module load Python/3.11.5-GCCcore-13.2.0 SciPy-bundle/2025.07-gfbf-2025b
  python3 10A-phylo_ssd.py --report REPORT --classify CLASSIFY --tree TREE \\
      --outdir FIGURES
"""

import argparse, csv, math, os, sys
from collections import Counter, defaultdict

try:
    import numpy as np
except ImportError:
    sys.exit("numpy required.  module load SciPy-bundle/2025.07-gfbf-2025b")
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
except ImportError:
    sys.exit("matplotlib required.\n"
             "  pip install --user matplotlib\n"
             "or module load a matplotlib-bearing bundle.")

RESIDUES = ["Tyr70", "Asp72", "Trp84", "Gly118", "Gly119", "Tyr121", "Glu199",
            "Ser200", "Ala201", "Trp279", "Phe288", "Phe290", "Glu327",
            "Phe330", "Tyr334", "His440"]
CANON = dict(zip(RESIDUES, "YDWGGYESAWFFEFYH"))


# ------------------------------------------------------------------ newick --

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
                sys.exit("newick parse error")
        j = i
        while i < len(text) and text[i] not in "(),:;": i += 1
        n.name = text[j:i].strip().strip("'\"")
        if i < len(text) and text[i] == ":":
            i += 1; j = i
            while i < len(text) and (text[i].isdigit() or text[i] in ".eE+-"): i += 1
            try: n.length = float(text[j:i])
            except ValueError: n.length = 0.0
        return n
    return node()


def prune(root, keep):
    """Keep only the named leaves, collapsing internal nodes and summing branch
    lengths so patristic distances are preserved."""
    def rec(n):
        if not n.children:
            return n if n.name in keep else None
        kids = [k for k in (rec(c) for c in n.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            kids[0].length += n.length
            return kids[0]
        m = Node(); m.name = n.name; m.length = n.length
        m.children = kids
        for k in kids: k.parent = m
        return m
    return rec(root)


def leaf_paths(root):
    out = {}
    def rec(n, path, dist):
        if not n.children:
            out[n.name] = (path, dist)
            return
        for c in n.children:
            rec(c, path + [n], dist + c.length)
    rec(root, [], 0.0)
    return out


def vcv_matrix(root, order):
    """Brownian-motion covariance: shared root-to-MRCA path length."""
    paths = leaf_paths(root)
    depth = {}
    def rec(n, d):
        depth[n] = d
        for c in n.children: rec(c, d + c.length)
    rec(root, 0.0)
    n = len(order)
    V = np.zeros((n, n))
    for i, a in enumerate(order):
        pa, da = paths[a]
        for j, b in enumerate(order):
            pb, _ = paths[b]
            common = [x for x in pa if x in pb]
            mrca = common[-1] if common else root
            V[i, j] = depth[mrca]
        V[i, i] = da
    return V


def lambda_ml(y, V):
    """ML estimate of Pagel's lambda under Brownian motion.

    lambda scales the off-diagonals of the covariance matrix: 0 means the trait
    is independent of the tree, 1 means pure Brownian evolution along it.
    """
    n = len(y)
    diag = np.diag(np.diag(V))
    off = V - diag

    def negll(lam):
        C = diag + lam * off
        try:
            L = np.linalg.cholesky(C + 1e-9 * np.eye(n))
        except np.linalg.LinAlgError:
            return 1e9
        Ci = np.linalg.inv(C + 1e-9 * np.eye(n))
        one = np.ones((n, 1))
        yv = y.reshape(-1, 1)
        # .item() rather than float(): numpy 2 refuses float() on a (1,1) array
        mu = ((one.T @ Ci @ yv) / (one.T @ Ci @ one)).item()
        r = yv - mu
        s2 = ((r.T @ Ci @ r) / n).item()
        if s2 <= 0:
            return 1e9
        logdet = 2 * np.sum(np.log(np.diag(L)))
        return 0.5 * (n * math.log(2 * math.pi * s2) + logdet + n)

    grid = np.linspace(0, 1, 201)
    lls = np.array([negll(l) for l in grid])
    best = grid[int(np.argmin(lls))]
    # local refinement
    lo, hi = max(0.0, best - 0.01), min(1.0, best + 0.01)
    fine = np.linspace(lo, hi, 81)
    ll2 = np.array([negll(l) for l in fine])
    lam = float(fine[int(np.argmin(ll2))])
    return lam, -negll(lam), -negll(0.0), -negll(1.0), grid, -lls


# ------------------------------------------------------------------- plots --

def draw_tree(ax, root, labels, colours, title):
    ypos, order = {}, []
    def assign(n):
        if not n.children:
            order.append(n); ypos[n] = len(order) - 1; return ypos[n]
        ys = [assign(c) for c in n.children]
        ypos[n] = sum(ys) / len(ys); return ypos[n]
    assign(root)
    xpos = {}
    def depth(n, d):
        xpos[n] = d
        for c in n.children: depth(c, d + c.length)
    depth(root, 0.0)
    for n in list(xpos):
        if n.parent is not None:
            ax.plot([xpos[n.parent], xpos[n]], [ypos[n], ypos[n]],
                    color="0.35", lw=0.8, zorder=1)
        if n.children:
            ys = [ypos[c] for c in n.children]
            ax.plot([xpos[n], xpos[n]], [min(ys), max(ys)],
                    color="0.35", lw=0.8, zorder=1)
    xm = max(xpos.values()) or 1.0
    for n in order:
        lab = labels.get(n.name, n.name)
        ax.scatter([xpos[n]], [ypos[n]], s=18, color=colours.get(n.name, "0.5"),
                   zorder=3, edgecolors="none")
        ax.text(xm * 1.02, ypos[n], lab, va="center", fontsize=6.5)
    ax.set_xlim(-0.01 * xm, xm * 1.65)
    ax.set_ylim(-1, len(order))
    ax.set_title(title, fontsize=10)
    ax.axis("off")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, help="REPORT/ with analysis_set.tsv")
    ap.add_argument("--classify", required=True)
    ap.add_argument("--tree", required=True, help="TREE/ or a Newick file")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--hc-level", type=float, default=5.0)
    ap.add_argument("--boot", type=int, default=2000)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    aset = list(csv.DictReader(
        open(os.path.join(a.report, "analysis_set.tsv"), newline=""),
        delimiter="\t"))
    calls = {r["sequence_id"]: r for r in csv.DictReader(
        open(os.path.join(a.classify, "active_site_calls.tsv"), newline=""),
        delimiter="\t")}
    tpath = (a.tree if os.path.isfile(a.tree)
             else os.path.join(a.tree, "ache.tree"))
    root = parse_newick(open(tpath).read())

    keep = {r["Sequence_id"] for r in aset}
    pr = prune(root, keep)
    if pr is None:
        sys.exit("no analysis sequences found in the tree")
    present = set(leaf_paths(pr))
    aset = [r for r in aset if r["Sequence_id"] in present]
    print(f"analysis set: {len(aset)} sequences, {len(present)} in the tree")

    classes = sorted({r["Class"] for r in aset})
    pal = plt.get_cmap("tab10")
    ccol = {c: pal(i % 10) for i, c in enumerate(classes)}
    labels = {r["Sequence_id"]:
              f"{r['Species_representative']}  ({float(r['log10_LC50']):.2f})"
              for r in aset}
    colours = {r["Sequence_id"]: ccol[r["Class"]] for r in aset}

    # ---- fig 1: tree ------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, max(5, 0.22 * len(aset))))
    draw_tree(ax, pr, labels, colours,
              "AChE gene tree, analysis set (tip label: log10 LC50 mg/L)")
    ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=ccol[c],
                                  label=c, markersize=5) for c in classes],
              loc="lower left", fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "fig1_tree.png"), dpi=200)
    fig.savefig(os.path.join(a.outdir, "fig1_tree.svg")); plt.close(fig)

    # ---- fig 2: active site ----------------------------------------------
    order = sorted(aset, key=lambda r: (r["Class"], r["Species_representative"]))
    M = np.zeros((len(order), len(RESIDUES)))
    txt = []
    for i, r in enumerate(order):
        c = calls.get(r["Sequence_id"], {})
        row = []
        for j, p in enumerate(RESIDUES):
            v = c.get(p, "-")
            row.append(v)
            M[i, j] = 0 if v == "-" else (1 if v == CANON[p] else 2)
        txt.append(row)
    fig, ax = plt.subplots(figsize=(8, max(4, 0.22 * len(order))))
    ax.imshow(M, aspect="auto", cmap=ListedColormap(["#dddddd", "#2c7fb8", "#e8703a"]),
              vmin=0, vmax=2)
    for i in range(len(order)):
        for j in range(len(RESIDUES)):
            ax.text(j, i, txt[i][j], ha="center", va="center", fontsize=5.5,
                    color="white" if M[i, j] else "0.4")
    ax.set_xticks(range(len(RESIDUES)))
    ax.set_xticklabels(RESIDUES, rotation=90, fontsize=7)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([r["Species_representative"] for r in order], fontsize=6)
    ax.set_title("Active-site residues, Torpedo numbering\n"
                 "blue = canonical, orange = substituted, grey = gap", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "fig2_active_site.png"), dpi=200)
    fig.savefig(os.path.join(a.outdir, "fig2_active_site.svg")); plt.close(fig)

    # ---- fig 3: domain completeness --------------------------------------
    gaps = [int(calls.get(r["Sequence_id"], {}).get("n_gaps_at_sites", 0) or 0)
            for r in order]
    arom = [int(calls.get(r["Sequence_id"], {}).get("n_gorge_aromatic", 0) or 0)
            for r in order]
    fig, axes = plt.subplots(1, 2, figsize=(9, max(4, 0.2 * len(order))))
    axes[0].barh(range(len(order)), [16 - g for g in gaps],
                 color=[ccol[r["Class"]] for r in order])
    axes[0].set_xlabel("diagnostic positions resolved (of 16)", fontsize=8)
    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels([r["Species_representative"] for r in order], fontsize=6)
    axes[0].invert_yaxis()
    axes[1].barh(range(len(order)), arom, color=[ccol[r["Class"]] for r in order])
    axes[1].axvline(8, ls="--", c="0.4", lw=0.8)
    axes[1].set_xlabel("aromatic residues lining the gorge", fontsize=8)
    axes[1].set_yticks([]); axes[1].invert_yaxis()
    fig.suptitle("Domain completeness and gorge aromaticity", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "fig3_coverage.png"), dpi=200)
    fig.savefig(os.path.join(a.outdir, "fig3_coverage.svg")); plt.close(fig)

    # ---- fig 4: Pagel's lambda -------------------------------------------
    ids = [r["Sequence_id"] for r in aset]
    y = np.array([float(r["log10_LC50"]) for r in aset])
    V = vcv_matrix(pr, ids)
    lam, ll, ll0, ll1, grid, curve = lambda_ml(y, V)
    lr0 = 2 * (ll - ll0)
    lr1 = 2 * (ll - ll1)
    print("\n== phylogenetic signal ==")
    print(f"  Pagel's lambda = {lam:.3f}")
    print(f"  logLik(lambda) {ll:.2f} | logLik(0) {ll0:.2f} | logLik(1) {ll1:.2f}")
    print(f"  LR vs lambda=0: {lr0:.2f}  (1 df; >3.84 rejects no signal)")
    print(f"  LR vs lambda=1: {lr1:.2f}  (1 df; >3.84 rejects Brownian)")
    fig, ax = plt.subplots(figsize=(5, 3.4))
    ax.plot(grid, curve, color="#2c7fb8")
    ax.axvline(lam, color="#e8703a", ls="--",
               label=f"ML $\\lambda$ = {lam:.3f}")
    ax.set_xlabel("Pagel's $\\lambda$"); ax.set_ylabel("log-likelihood")
    ax.set_title("Phylogenetic signal in log10 LC50", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "fig4_lambda.png"), dpi=200)
    fig.savefig(os.path.join(a.outdir, "fig4_lambda.svg")); plt.close(fig)

    # ---- fig 5: SSD -------------------------------------------------------
    ys = np.sort(y)
    n = len(ys)
    emp = (np.arange(1, n + 1) - 0.5) / n
    mu, sd = float(np.mean(ys)), float(np.std(ys, ddof=1))

    def norm_ppf(p):
        # Acklam rational approximation; adequate for HC5 reporting
        a_ = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
              1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
        b_ = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
              6.680131188771972e+01, -1.328068155288572e+01]
        c_ = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
              -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
        d_ = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
              3.754408661907416e+00]
        pl = 0.02425
        if p < pl:
            q = math.sqrt(-2 * math.log(p))
            return (((((c_[0]*q+c_[1])*q+c_[2])*q+c_[3])*q+c_[4])*q+c_[5]) / \
                   ((((d_[0]*q+d_[1])*q+d_[2])*q+d_[3])*q+1)
        if p > 1 - pl:
            q = math.sqrt(-2 * math.log(1 - p))
            return -(((((c_[0]*q+c_[1])*q+c_[2])*q+c_[3])*q+c_[4])*q+c_[5]) / \
                    ((((d_[0]*q+d_[1])*q+d_[2])*q+d_[3])*q+1)
        q = p - 0.5; r = q * q
        return (((((a_[0]*r+a_[1])*r+a_[2])*r+a_[3])*r+a_[4])*r+a_[5])*q / \
               (((((b_[0]*r+b_[1])*r+b_[2])*r+b_[3])*r+b_[4])*r+1)

    z5 = norm_ppf(a.hc_level / 100.0)
    hc5_log = mu + z5 * sd
    rng = np.random.default_rng(0)
    boots = []
    for _ in range(a.boot):
        s = rng.choice(ys, size=n, replace=True)
        boots.append(np.mean(s) + z5 * np.std(s, ddof=1))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    print(f"\n== SSD (log-normal, n={n}) ==")
    print(f"  mean {mu:.3f}, sd {sd:.3f} (log10 mg/L)")
    print(f"  HC{a.hc_level:g} = {10**hc5_log:.4g} mg/L "
          f"(95% CI {10**lo:.4g} to {10**hi:.4g})")

    xs = np.linspace(ys.min() - 0.7, ys.max() + 0.7, 400)
    cdf = np.array([0.5 * (1 + math.erf((x - mu) / (sd * math.sqrt(2))))
                    for x in xs])
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.plot(xs, cdf * 100, color="0.25", lw=1.6, zorder=2,
            label="log-normal fit")
    # Points carry class as colour and order as marker, so the SSD shows WHICH
    # taxa occupy the sensitive tail rather than only that a tail exists. An
    # SSD without taxonomy is not interpretable for risk assessment: the whole
    # question is which groups sit below the HC5.
    orders = sorted({r.get("Order", "") for r in aset if r.get("Order")})
    marks = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "*", "h", "p"]
    omark = {o: marks[i % len(marks)] for i, o in enumerate(orders)}
    for r in aset:
        v = float(r["log10_LC50"])
        pct = 100 * (np.searchsorted(ys, v) - 0.5) / n
        ax.scatter([v], [pct], color=ccol[r["Class"]], s=46, zorder=4,
                   marker=omark.get(r.get("Order", ""), "o"),
                   edgecolors="0.25", linewidths=0.4)
    # label the five most sensitive species; those drive the HC5
    for r in sorted(aset, key=lambda x: float(x["log10_LC50"]))[:5]:
        v = float(r["log10_LC50"])
        pct = 100 * (np.searchsorted(ys, v) - 0.5) / n
        ax.annotate(r["Species_representative"], (v, pct),
                    textcoords="offset points", xytext=(7, -2), fontsize=6.5,
                    style="italic", color="0.3")
    ax.axhline(a.hc_level, color="0.5", ls=":", lw=0.9)
    ax.axvline(hc5_log, color="#c0392b", ls="--", lw=1.3,
               label=f"HC{a.hc_level:g} = {10**hc5_log:.3g} mg/L")
    ax.fill_betweenx([0, 100], lo, hi, color="#c0392b", alpha=0.10,
                     label="95% CI (bootstrap)", zorder=1)
    ax.set_xlabel("log10 LC50 (mg/L)")
    ax.set_ylabel("species affected (%)")
    ax.set_title(f"Species sensitivity distribution, chlorpyrifos (n={n})",
                 fontsize=10)
    ax.set_ylim(0, 100)
    h1 = [plt.Line2D([], [], marker="o", ls="", color=ccol[c], label=c,
                     markersize=6, markeredgecolor="0.25") for c in classes]
    h2 = [plt.Line2D([], [], marker=omark[o], ls="", color="0.55", label=o,
                     markersize=5.5, markeredgecolor="0.25") for o in orders]
    leg1 = ax.legend(handles=h1, title="class", fontsize=7, title_fontsize=7,
                     loc="upper left", frameon=False)
    ax.add_artist(leg1)
    if orders:
        ax.legend(handles=h2, title="order", fontsize=6, title_fontsize=6.5,
                  loc="lower right", frameon=False, ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "fig5_ssd.png"), dpi=200)
    fig.savefig(os.path.join(a.outdir, "fig5_ssd.svg")); plt.close(fig)

    # per-class contribution to the sensitive tail
    below = [r for r in aset if float(r["log10_LC50"]) <= hc5_log]
    print(f"  species at or below the HC{a.hc_level:g}: {len(below)}")
    for c, k in Counter(r["Class"] for r in below).most_common():
        print(f"    {c}: {k}")

    # ---- fig 6: by class and clade ---------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, key, ttl in ((axes[0], "Class", "by taxonomic class"),
                         (axes[1], "AChE_clade", "by AChE clade")):
        groups = sorted({r[key] for r in aset})
        data = [[float(r["log10_LC50"]) for r in aset if r[key] == g]
                for g in groups]
        ax.boxplot(data, vert=True, widths=0.5,
                   medianprops=dict(color="#e8703a"))
        for i, d in enumerate(data, 1):
            ax.scatter(np.full(len(d), i) + rng.normal(0, 0.05, len(d)), d,
                       s=16, alpha=0.75, color="#2c7fb8", edgecolors="none")
        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels([f"{g}\n(n={len(d)})" for g, d in zip(groups, data)],
                           fontsize=7, rotation=30, ha="right")
        ax.set_ylabel("log10 LC50 (mg/L)"); ax.set_title(ttl, fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "fig6_by_class.png"), dpi=200)
    fig.savefig(os.path.join(a.outdir, "fig6_by_class.svg")); plt.close(fig)

    with open(os.path.join(a.outdir, "stats.txt"), "w") as fh:
        fh.write(f"n\t{len(aset)}\nlambda\t{lam:.4f}\nlogLik\t{ll:.4f}\n"
                 f"logLik_lambda0\t{ll0:.4f}\nlogLik_lambda1\t{ll1:.4f}\n"
                 f"LR_vs_0\t{lr0:.4f}\nLR_vs_1\t{lr1:.4f}\n"
                 f"ssd_mean_log10\t{mu:.4f}\nssd_sd_log10\t{sd:.4f}\n"
                 f"HC{a.hc_level:g}_mg_L\t{10**hc5_log:.6g}\n"
                 f"HC_lo95\t{10**lo:.6g}\nHC_hi95\t{10**hi:.6g}\n")

    print(f"\nfigures and stats.txt -> {a.outdir}/")
    if lam > 0.8:
        print(f"\nlambda = {lam:.2f}. Most of the interspecies variance is "
              "phylogenetic. Fit the structural predictors as PGLS on the "
              "residual and report the effect size honestly; a weak effect "
              "here is a real finding, not a failed analysis.")
    elif lam < 0.2:
        print(f"\nlambda = {lam:.2f}. Little phylogenetic signal, so ordinary "
              "regression is defensible and there is real residual variance "
              "for a structural predictor to explain.")


if __name__ == "__main__":
    main()
