"""netZooPy reference networks for the PANDA tests (shared with accNetZoo).

Run with the BEAT virtualenv (netZooPy installed):
  python make_reference.py <output folder>
For each case it writes the inputs and netZooPy's outputs as CSV with 17
significant digits, so R reads the exact doubles.
"""
import json, os, sys, tempfile
import numpy as np
import pandas as pd
from netZooPy.panda import Panda
from netZooPy.lioness import Lioness
from netZooPy.panda import calculations as calc

out_root = sys.argv[1]
FMT = "%.17g"


def panda_all(corr, ppi, motif, alpha, threshold=0.001):
    """compute_panda_cpu, returning all three networks and the step count."""
    motif, ppi, corr = motif.copy(), ppi.copy(), corr.copy()
    nt, ng = motif.shape
    step, hamming = 0, 1
    while hamming > threshold:
        W = 0.5 * (calc.t_function(ppi, motif) + calc.t_function(motif, corr))
        hamming = np.abs(motif - W).mean()
        motif *= 1 - alpha
        motif += alpha * W
        p = calc.update_diagonal(calc.t_function(motif), nt, alpha, step)
        ppi *= 1 - alpha
        ppi += alpha * p
        c = calc.update_diagonal(calc.t_function(motif.T), ng, alpha, step)
        corr *= 1 - alpha
        corr += alpha * c
        step += 1
    return motif, ppi, corr, step


def save(df, path, **kw):
    df.to_csv(path, float_format=FMT, **kw)


def case(name, expr, motif, ppi, mode="intersection", alpha=0.1, lioness=(), full=True):
    d = os.path.join(out_root, name)
    os.makedirs(d, exist_ok=True)
    if expr is not None:
        save(expr, os.path.join(d, "expr.csv"))
    save(motif, os.path.join(d, "motif.csv"), index=False)
    if ppi is not None:
        save(ppi, os.path.join(d, "ppi.csv"), index=False)
    m = motif.copy(); m.columns = [0, 1, 2]
    p = None
    if ppi is not None:
        p = ppi.copy(); p.columns = [0, 1, 2]
    pa = Panda(expr, m, p, computing="cpu", precision="double", save_memory=False,
               keep_expression_matrix=len(lioness) > 0, modeProcess=mode, alpha=alpha)
    reg = pa.panda_network
    save(reg, os.path.join(d, "regulatory.csv"))
    motif2, ppi2, corr2, steps = panda_all(pa.correlation_matrix, pa.ppi_matrix,
                                           pa.motif_matrix, alpha)
    assert np.array_equal(motif2, reg.to_numpy()), "copy of compute_panda differs"
    if full:
        save(pd.DataFrame(ppi2, index=reg.index, columns=reg.index),
             os.path.join(d, "cooperative.csv"))
        save(pd.DataFrame(corr2, index=reg.columns, columns=reg.columns),
             os.path.join(d, "coexpression.csv"))
    meta = {"mode": mode, "alpha": alpha, "iterations": steps, "lioness": list(lioness)}
    if lioness:
        li = Lioness(pa, computing="cpu", precision="double", ncores=1,
                     subset_numbers=list(lioness), save_dir=tempfile.mkdtemp(),
                     alpha=alpha)
        res = li.export_lioness_results
        res.columns = [str(c) for c in res.columns]
        save(res, os.path.join(d, "lioness.csv"), index=False)
    with open(os.path.join(d, "meta.json"), "w") as f:
        json.dump(meta, f)


rng = np.random.default_rng(1)


def make(nt, ng, ns, density=0.15, extra_genes=0, seed=None):
    tfs = [f"TF{i:02d}" for i in range(nt)]
    genes = [f"G{i:04d}" for i in range(ng)]
    expr = pd.DataFrame(rng.normal(size=(ng, ns)) + rng.normal(size=(ng, 1)) * 0.5,
                        index=genes, columns=[f"S{i:02d}" for i in range(ns)])
    # some structure: TF-driven modules
    act = rng.normal(size=(nt, ns))
    for t in range(nt):
        tg = rng.choice(ng, size=max(3, ng // 20), replace=False)
        expr.iloc[tg, :] += 0.8 * act[t]
    rows = []
    for t in tfs:
        for g in rng.choice(genes + [f"X{i}" for i in range(extra_genes)],
                            size=int(density * ng), replace=False):
            rows.append((t, g, float(rng.choice([1.0, 0.5, rng.random()]))))
    motif = pd.DataFrame(rows, columns=["source", "target", "weight"])
    pairs = [(tfs[i], tfs[j], float(rng.random())) for i in range(nt) for j in range(i + 1, nt)
             if rng.random() < 0.15]
    ppi = pd.DataFrame(pairs, columns=["a", "b", "w"])
    return expr, motif, ppi, tfs, genes


# 1: plain intersection mode, weighted motif with repeated pairs, LIONESS
expr, motif, ppi, tfs, genes = make(25, 120, 20)
motif = pd.concat([motif, motif.iloc[:15].assign(weight=0.25)], ignore_index=True)
case("basic", expr, motif, ppi, lioness=(0, 7, 19))

# 2: tricky inputs in union mode: genes without expression (zero rows), a
# constant gene, a TF with no targets (only in the PPI), a pair given in both
# directions with different weights, self-pairs, unsorted names
expr, motif, ppi, tfs, genes = make(15, 80, 12, extra_genes=10)
expr.iloc[3, :] = 2.0
ppi = pd.concat([ppi, pd.DataFrame([(tfs[1], tfs[0], 0.9), (tfs[0], tfs[1], 0.3),
                                     (tfs[2], tfs[2], 0.7), ("TFonlyPPI", tfs[3], 1.0)],
                                    columns=["a", "b", "w"])], ignore_index=True)
expr = expr.iloc[rng.permutation(len(expr))]
case("union", expr, motif, ppi, mode="union", alpha=0.2, lioness=(0, 11))
case("union_intersect", expr, motif, ppi, mode="intersection", lioness=(5,))

# 3: no PPI (identity) and no expression (identity)
expr, motif, ppi, tfs, genes = make(10, 50, 8)
case("noppi", expr, motif, None, lioness=(2,))
case("noexpr", None, motif, ppi)

# 4: larger, to cross block boundaries; only the regulatory network and two
# LIONESS networks are stored
expr, motif, ppi, tfs, genes = make(60, 1500, 40, density=0.05)
case("large", expr, motif, ppi, lioness=(0, 39), full=False)
print("done")
