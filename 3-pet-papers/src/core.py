"""Core numerics for the statistical review of Goold et al. (2016), R. Soc. Open Sci. 3:160268.

Re-derives the authors' pipeline from their own Rscript.R:
  polychoric correlations (cor_auto)  ->  EBICglasso(gamma=0.65)  ->  centrality
Implemented here because scikit-learn / networkx / an R bridge are not available in
this sandbox; every routine follows the conventional algorithm and is checked against
the paper's own reported quantities (edge counts 55 / 70).
"""
import numpy as np
from scipy import stats, optimize

# ---------------------------------------------------------------- bivariate normal CDF
_GLX, _GLW = np.polynomial.legendre.leggauss(40)


def _bvn(h, k, rho):
    """P(X<=h, Y<=k) for standard bivariate normal with correlation rho.

    Plackett / Sheppard identity: Phi2(h,k,r) = Phi(h)Phi(k) + int_0^r phi2(h,k,t) dt,
    evaluated by 40-point Gauss-Legendre quadrature.  h, k broadcast; rho scalar.
    """
    h = np.asarray(h, float)
    k = np.asarray(k, float)
    out = stats.norm.cdf(h) * stats.norm.cdf(k)
    if rho == 0.0:
        return out
    t = 0.5 * rho * (_GLX + 1.0)                      # nodes on [0, rho]
    w = 0.5 * rho * _GLW
    hh = h[..., None]
    kk = k[..., None]
    om = 1.0 - t ** 2
    dens = np.exp(-(hh ** 2 - 2 * t * hh * kk + kk ** 2) / (2 * om)) / (2 * np.pi * np.sqrt(om))
    dens = np.where(np.isfinite(dens), dens, 0.0)
    return out + (dens * w).sum(axis=-1)


def _thresholds(counts):
    """Normal thresholds from marginal counts; returns vector of length ncat+1 (-inf..inf)."""
    p = np.cumsum(counts) / counts.sum()
    tau = stats.norm.ppf(np.clip(p[:-1], 1e-12, 1 - 1e-12))
    return np.concatenate(([-np.inf], tau, [np.inf]))


def polychoric(x, y):
    """Two-step maximum-likelihood polychoric correlation for two ordinal vectors."""
    x = np.asarray(x, int)
    y = np.asarray(y, int)
    cx = np.unique(x)
    cy = np.unique(y)
    if len(cx) < 2 or len(cy) < 2:
        return 0.0
    xi = np.searchsorted(cx, x)
    yi = np.searchsorted(cy, y)
    tab = np.zeros((len(cx), len(cy)))
    np.add.at(tab, (xi, yi), 1.0)
    a = _thresholds(tab.sum(axis=1))
    b = _thresholds(tab.sum(axis=0))
    A, B = np.meshgrid(a, b, indexing="ij")

    def nll(rho):
        P = _bvn(A, B, float(rho))
        cell = P[1:, 1:] - P[:-1, 1:] - P[1:, :-1] + P[:-1, :-1]
        cell = np.clip(cell, 1e-12, None)
        return -(tab * np.log(cell)).sum()

    r = optimize.minimize_scalar(nll, bounds=(-0.999, 0.999), method="bounded",
                                 options={"xatol": 1e-5}).x
    return float(r)


def polychoric_matrix(X):
    """Polychoric correlation matrix of an (n x p) integer ordinal array."""
    X = np.asarray(X, int)
    p = X.shape[1]
    R = np.eye(p)
    for i in range(p):
        for j in range(i + 1, p):
            R[i, j] = R[j, i] = polychoric(X[:, i], X[:, j])
    return _force_pd(R)


def _force_pd(R, eig_tol=1e-6, posd_tol=1e-8, max_iter=200, conv_tol=1e-7):
    """Nearest positive-definite correlation matrix, Higham (2002) alternating projections.

    Mirrors Matrix::nearPD(corr = TRUE), which is what qgraph::cor_auto(forcePD = TRUE)
    calls: project onto the PSD cone (dropping eigenvalues <= eig_tol * lambda_max),
    then onto unit diagonal, with Dykstra correction; finish by flooring the spectrum
    at posd_tol * lambda_max.
    """
    R = np.asarray(R, float)
    if np.linalg.eigvalsh(R).min() > posd_tol:
        return R
    n = R.shape[0]
    X = R.copy()
    D_S = np.zeros_like(X)
    for _ in range(max_iter):
        Y = X.copy()
        Rk = Y - D_S
        w, V = np.linalg.eigh(Rk)
        keep = w > eig_tol * w.max()
        X = (V[:, keep] * w[keep]) @ V[:, keep].T
        D_S = X - Rk
        np.fill_diagonal(X, 1.0)
        if np.abs(Y - X).max() / max(np.abs(Y).max(), 1e-12) < conv_tol:
            break
    w, V = np.linalg.eigh(X)
    eps = posd_tol * abs(w.max())
    if w.min() < eps:
        w = np.clip(w, eps, None)
        X = (V * w) @ V.T
    d = np.sqrt(np.diag(X))
    X = X / np.outer(d, d)
    return (X + X.T) / 2


# ---------------------------------------------------------------- graphical lasso
def glasso(S, rho, tol=1e-5, max_iter=100, warm=None):
    """Friedman/Hastie/Tibshirani (2008) block coordinate-descent graphical lasso.

    Minimises -log det(Theta) + tr(S Theta) + rho * ||Theta||_1 over off-diagonals
    (penalize.diagonal = FALSE, as called by qgraph::EBICglasso).
    """
    p = S.shape[0]
    if warm is None:
        W = S.copy()
        B = np.zeros((p, p))
    else:
        W, B = warm[0].copy(), warm[1].copy()
    idx = np.arange(p)
    thresh = tol * np.abs(S[~np.eye(p, dtype=bool)]).mean()
    for _ in range(max_iter):
        W_old = W.copy()
        for j in range(p):
            m = idx != j
            W11 = W[np.ix_(m, m)]
            s12 = S[m, j]
            beta = B[m, j].copy()
            d = np.diag(W11).copy()
            d[d < 1e-12] = 1e-12
            r = s12 - W11 @ beta                      # residual, updated in place
            for _inner in range(200):
                delta = 0.0
                for k in range(p - 1):
                    bk = beta[k]
                    rk = r[k] + d[k] * bk
                    new = np.sign(rk) * max(abs(rk) - rho, 0.0) / d[k]
                    if new != bk:
                        r -= W11[:, k] * (new - bk)
                        beta[k] = new
                        delta = max(delta, abs(new - bk))
                if delta < 1e-7:
                    break
            B[m, j] = beta
            w12 = W11 @ beta
            W[m, j] = w12
            W[j, m] = w12
        if np.abs(W - W_old).mean() < thresh:
            break
    Theta = np.zeros((p, p))
    for j in range(p):
        m = idx != j
        beta = B[m, j]
        tjj = 1.0 / max(W[j, j] - W[m, j] @ beta, 1e-12)
        Theta[j, j] = tjj
        Theta[m, j] = -beta * tjj
    return (Theta + Theta.T) / 2


def wi2net(Theta):
    """Precision matrix -> partial correlation network (qgraph::wi2net)."""
    d = np.sqrt(np.abs(np.diag(Theta)))
    P = -Theta / np.outer(d, d)
    np.fill_diagonal(P, 0.0)
    P[np.abs(P) < 1e-10] = 0.0
    return P


def _W_from(S, Theta):
    return np.linalg.inv(Theta)


def _B_from(Theta):
    """Recover the lasso coefficients beta_j = -Theta[-j,j]/Theta[jj] for warm starting."""
    p = Theta.shape[0]
    B = np.zeros((p, p))
    for j in range(p):
        m = np.arange(p) != j
        B[m, j] = -Theta[m, j] / Theta[j, j]
    return B


def ebic_glasso(S, n, gamma=0.65, nlambda=100, lambda_min_ratio=0.01, return_path=False):
    """qgraph::EBICglasso -- glasso over a lambda path, lambda chosen by minimum EBIC."""
    p = S.shape[0]
    off = S - np.eye(p)
    lam_max = max(off.max(), -off.min())
    lams = np.exp(np.linspace(np.log(lam_max * lambda_min_ratio), np.log(lam_max), nlambda))
    best = None
    path = []
    warm = None
    for lam in lams[::-1]:                 # descend the path with warm starts
        Theta = glasso(S, lam, warm=warm)
        warm = (_W_from(S, Theta), _B_from(Theta))
        net = wi2net(Theta)
        E = int((np.abs(np.triu(net, 1)) > 0).sum())
        sign_, logdet = np.linalg.slogdet(Theta)
        if sign_ <= 0:
            continue
        L = n / 2.0 * (logdet - np.trace(S @ Theta))
        ebic = -2 * L + E * np.log(n) + 4 * E * gamma * np.log(p)
        path.append((lam, E, ebic))
        if best is None or ebic < best["ebic"]:
            best = {"lambda": lam, "ebic": ebic, "net": net, "Theta": Theta, "n_edges": E}
    if return_path:
        best["path"] = np.array(path)
    return best


# ---------------------------------------------------------------- centrality
def strength(net):
    """qgraph strength centrality: sum of absolute edge weights per node."""
    return np.abs(net).sum(axis=1)


def betweenness(net):
    """Weighted betweenness (Brandes 2001) on edge lengths 1/|w|, undirected.

    Matches the igraph/qgraph convention: each unordered pair counted once.
    """
    import heapq
    p = net.shape[0]
    A = np.abs(net)
    with np.errstate(divide="ignore"):
        length = np.where(A > 0, 1.0 / np.where(A > 0, A, 1.0), np.inf)
    nbrs = [np.flatnonzero(A[i] > 0) for i in range(p)]
    CB = np.zeros(p)
    for s in range(p):
        S_stack, P = [], [[] for _ in range(p)]
        sigma = np.zeros(p)
        sigma[s] = 1.0
        seen = {s: 0.0}
        Q = [(0.0, s)]
        done = np.zeros(p, bool)
        while Q:
            d, v = heapq.heappop(Q)
            if done[v]:
                continue
            done[v] = True
            S_stack.append(v)
            for w in nbrs[v]:
                vw = d + length[v, w]
                if done[w]:
                    continue
                if w not in seen or vw < seen[w] - 1e-12:
                    seen[w] = vw
                    heapq.heappush(Q, (vw, w))
                    sigma[w] = sigma[v]
                    P[w] = [v]
                elif abs(vw - seen[w]) < 1e-12:
                    sigma[w] += sigma[v]
                    P[w].append(v)
        delta = np.zeros(p)
        for w in reversed(S_stack):
            for v in P[w]:
                delta[v] += sigma[v] / sigma[w] * (1 + delta[w])
            if w != s:
                CB[w] += delta[w]
    return CB / 2.0


def network_from_data(X, n=None, gamma=0.65):
    """Full pipeline: ordinal data -> polychoric -> EBICglasso -> network + centrality."""
    X = np.asarray(X, int)
    n = X.shape[0] if n is None else n
    R = polychoric_matrix(X)
    fit = ebic_glasso(R, n=n, gamma=gamma)
    fit["R"] = R
    fit["strength"] = strength(fit["net"])
    fit["betweenness"] = betweenness(fit["net"])
    return fit


# ---------------------------------------------------------------- misc
def cliffs_delta(a, b):
    """Cliff's delta: P(a>b) - P(a<b), computed by rank identity."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    na, nb = len(a), len(b)
    allv = np.concatenate([a, b])
    r = stats.rankdata(allv)
    ra = r[:na].sum()
    u = ra - na * (na + 1) / 2.0          # Mann-Whitney U for a
    return float(2 * u / (na * nb) - 1)
