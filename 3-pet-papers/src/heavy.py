"""Resampling workload for the Goold et al. (2016) review.  Run as a script.

Writes one .npz per job into goold_2016/results/ so a kernel restart loses nothing.
Seeds are fixed; every job records its own seed and replicate count.
"""
import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import core  # noqa: E402
import prep  # noqa: E402

OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

NLAMBDA_RESAMPLE = 25      # coarser than the 100 used for the two reported networks
LAM_MIN_RATIO_RESAMPLE = 0.05   # reduced resampling budget; see report section 10
GAMMA = 0.65
IU = np.triu_indices(20, 1)

SHARD, NSHARDS = 0, 1          # set from the command line


def _shard_map(fn, tasks):
    """Run only this process's slice of `tasks` (multiprocessing pipes are blocked here)."""
    return [fn(t) for i, t in enumerate(tasks) if i % NSHARDS == SHARD]


def _out(name):
    return os.path.join(OUT, f"{name}_s{SHARD}of{NSHARDS}.npz")


def fit(X, nlambda=NLAMBDA_RESAMPLE):
    R = core.polychoric_matrix(X)
    f = core.ebic_glasso(R, n=X.shape[0], gamma=GAMMA, nlambda=nlambda,
                         lambda_min_ratio=LAM_MIN_RATIO_RESAMPLE)
    net = f["net"]
    return net, f["n_edges"]


# ----------------------------------------------------------------- jobs
def _boot_one(args):
    X, seed = args
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, X.shape[0], X.shape[0])
    try:
        net, E = fit(X[idx])
    except Exception:
        return None
    return (E, core.strength(net), core.betweenness(net), net[IU])


def job_bootstrap(name, X, n_boot, seed0):
    tasks = [(X, seed0 + i) for i in range(n_boot)]
    res = [r for r in _shard_map(_boot_one, tasks) if r is not None]
    np.savez_compressed(
        _out(f"boot_{name}"),
        n_edges=np.array([r[0] for r in res]),
        strength=np.array([r[1] for r in res]),
        betweenness=np.array([r[2] for r in res]),
        edges=np.array([r[3] for r in res]),
        n_boot=len(res), seed0=seed0, nlambda=NLAMBDA_RESAMPLE,
    )
    return len(res)


def _case_one(args):
    X, frac, seed = args
    rng = np.random.default_rng(seed)
    m = max(10, int(round(frac * X.shape[0])))
    idx = rng.choice(X.shape[0], m, replace=False)
    try:
        net, E = fit(X[idx])
    except Exception:
        return None
    return (frac, m, core.strength(net), core.betweenness(net))


def job_casedrop(name, X, fracs, reps, seed0):
    tasks = [(X, f, seed0 + 1000 * k + i) for k, f in enumerate(fracs) for i in range(reps)]
    res = [r for r in _shard_map(_case_one, tasks) if r is not None]
    np.savez_compressed(
        _out(f"casedrop_{name}"),
        frac=np.array([r[0] for r in res]), m=np.array([r[1] for r in res]),
        strength=np.array([r[2] for r in res]), betweenness=np.array([r[3] for r in res]),
        reps=reps, seed0=seed0,
    )
    return len(res)


def _node_one(args):
    X, k, seed = args
    rng = np.random.default_rng(seed)
    keep = np.sort(rng.choice(X.shape[1], k, replace=False))
    Xi = X[:, keep]
    try:
        R = core.polychoric_matrix(Xi)
        f = core.ebic_glasso(R, n=Xi.shape[0], gamma=GAMMA, nlambda=NLAMBDA_RESAMPLE,
                             lambda_min_ratio=LAM_MIN_RATIO_RESAMPLE)
        net = f["net"]
    except Exception:
        return None
    s = np.full(X.shape[1], np.nan)
    b = np.full(X.shape[1], np.nan)
    s[keep] = core.strength(net)
    b[keep] = core.betweenness(net)
    return (k, s, b)


def job_nodedrop(name, X, ks, reps, seed0):
    tasks = [(X, k, seed0 + 1000 * j + i) for j, k in enumerate(ks) for i in range(reps)]
    res = [r for r in _shard_map(_node_one, tasks) if r is not None]
    np.savez_compressed(
        _out(f"nodedrop_{name}"),
        k=np.array([r[0] for r in res]),
        strength=np.array([r[1] for r in res]), betweenness=np.array([r[2] for r in res]),
        reps=reps, seed0=seed0,
    )
    return len(res)


def _nct_one(args):
    X, n1, seed = args
    rng = np.random.default_rng(seed)
    perm = rng.permutation(X.shape[0])
    A, B = X[perm[:n1]], X[perm[n1:]]
    try:
        na, _ = fit(A)
        nb, _ = fit(B)
    except Exception:
        return None
    return (na[IU] - nb[IU],
            core.strength(na) - core.strength(nb),
            core.betweenness(na) - core.betweenness(nb))


def job_nct(X, n1, n_perm, seed0):
    """Permutation test of patrol-vs-detection network difference (van Borkulo-style NCT)."""
    tasks = [(X, n1, seed0 + i) for i in range(n_perm)]
    res = [r for r in _shard_map(_nct_one, tasks) if r is not None]
    np.savez_compressed(
        _out("nct"),
        d_edge=np.array([r[0] for r in res]),
        d_strength=np.array([r[1] for r in res]),
        d_betweenness=np.array([r[2] for r in res]),
        n_perm=len(res), n1=n1, seed0=seed0,
    )
    return len(res)


def _maxstat_one(args):
    X, seed = args
    rng = np.random.default_rng(seed)
    Xp = np.column_stack([rng.permutation(X[:, j]) for j in range(X.shape[1])])
    try:
        R = core.polychoric_matrix(Xp)
        P = core.wi2net(np.linalg.inv(R))
    except Exception:
        return None
    return np.abs(R[IU]), np.abs(P[IU])


def job_maxstat(name, X, n_perm, seed0):
    """Null distribution of |partial correlation| over all 190 node pairs, columns
    permuted independently so every pairwise association is null by construction."""
    tasks = [(X, seed0 + i) for i in range(n_perm)]
    res = [r for r in _shard_map(_maxstat_one, tasks) if r is not None]
    np.savez_compressed(_out(f"maxstat_{name}"),
                        null_zero_order=np.array([r[0] for r in res]),
                        null_partial=np.array([r[1] for r in res]),
                        n_perm=len(res), seed0=seed0)
    return len(res)


if __name__ == "__main__":
    SHARD, NSHARDS = int(sys.argv[1]), int(sys.argv[2])
    raw, mi, _, _ = prep.load_all()
    patrol, detect = prep.network_inputs(mi)
    X_all = np.vstack([patrol, detect])
    FRACS = [0.25, 0.33, 0.41, 0.49, 0.57, 0.65, 0.73, 0.81, 0.89, 0.95]
    KS = list(range(3, 20, 2))
    t0 = time.time()
    log = []

    def run(label, fn, *a):
        t = time.time()
        n = fn(*a)
        log.append(f"{label}: {n} reps in {time.time()-t:.0f}s (total {time.time()-t0:.0f}s)")
        print(log[-1], flush=True)

    run("bootstrap patrol", job_bootstrap, "patrol", patrol, 100, 10001)
    run("bootstrap detect", job_bootstrap, "detect", detect, 100, 20001)
    run("nct", job_nct, X_all, 116, 100, 70001)
    run("maxstat patrol", job_maxstat, "patrol", patrol, 300, 90001)
    run("maxstat detect", job_maxstat, "detect", detect, 300, 91001)
    run("casedrop patrol", job_casedrop, "patrol", patrol, FRACS, 6, 30001)
    run("casedrop detect", job_casedrop, "detect", detect, FRACS, 6, 40001)
    run("nodedrop patrol", job_nodedrop, "patrol", patrol, KS, 6, 50001)
    run("nodedrop detect", job_nodedrop, "detect", detect, KS, 6, 60001)
    with open(os.path.join(OUT, f"heavy_log_s{SHARD}.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
