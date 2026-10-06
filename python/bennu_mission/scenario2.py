"""Scenario 2: direct two-impulse heliocentric transfer Earth -> target.

The transfer model is vectorized: evaluate_transfers() handles any number of
candidate solutions at once, so the grid search and the population of the
evolutionary optimizer are evaluated without Python loops.
"""

from types import SimpleNamespace as NS

import numpy as np
from scipy.optimize import differential_evolution, minimize

from .core import TWO_PI, _perifocal_axes, kep2car, time_of_flight


def evaluate_transfers(X, cfg):
    """Delta-v and time of flight for an array X (N, 3) of [nu_dep, nu_arr, argp_t].

    Returns a dict of arrays (N,) (vectors: (N, 3)). Same model and penalties
    as evaluateTransfer.m: 'valid' is False, with a positive 'penalty', when
    the transfer is not an ellipse or its perihelion is too low.
    """
    s = cfg.sc2
    mu = cfg.const.mu_sun
    X = np.mod(np.atleast_2d(np.asarray(X, float)), TWO_PI)
    n = X.shape[0]
    oE, oT = s.earth, s.target.orbit
    r1, v_earth = kep2car(oE.a, oE.e, oE.i, oE.raan, oE.argp, X[:, 0], mu)
    r2, v_target = kep2car(oT.a, oT.e, oT.i, oT.raan, oT.argp, X[:, 1], mu)
    r1n = np.linalg.norm(r1, axis=1)
    r2n = np.linalg.norm(r2, axis=1)

    out = {k: np.full(n, np.nan) for k in ("dv", "dv1", "dv2", "tof_days", "a_t", "e_t", "i_t",
                                           "raan_t", "nu1_t", "nu2_t", "arc_deg", "rp_t")}
    out.update(x=X, argp_t=X[:, 2], r1=r1, r2=r2, v_earth=v_earth, v_target=v_target,
               v_t1=np.full((n, 3), np.nan), v_t2=np.full((n, 3), np.nan),
               J=np.full(n, np.inf), valid=np.zeros(n, bool), penalty=np.zeros(n))

    with np.errstate(divide="ignore", invalid="ignore"):
        # transfer plane: contains r1 and r2 (prograde when long arcs are allowed)
        h = np.cross(r1, r2)
        hn = np.linalg.norm(h, axis=1)
        ok = hn >= 1e-8 * r1n * r2n
        h = h / hn[:, None]
        if s.allow_long_arc:
            h[h[:, 2] < 0] *= -1
        i_t = np.arccos(np.clip(h[:, 2], -1, 1))
        node = np.stack([-h[:, 1], h[:, 0], np.zeros(n)], 1)        # z x h
        nn = np.linalg.norm(node, axis=1)
        raan_t = np.where(nn < 1e-12, 0.0, np.mod(np.arctan2(node[:, 1] / nn, node[:, 0] / nn), TWO_PI))

        # true anomalies of r1 and r2 on the transfer orbit (perifocal components)
        p_axis, q_axis, _ = _perifocal_axes(raan_t, i_t, X[:, 2])
        nu1 = np.mod(np.arctan2(np.sum(q_axis * r1, 1), np.sum(p_axis * r1, 1)), TWO_PI)
        nu2 = np.mod(np.arctan2(np.sum(q_axis * r2, 1), np.sum(p_axis * r2, 1)), TWO_PI)

        # conic r = p / (1 + e cos(nu)) written at r1 and r2
        den = r1n * np.cos(nu1) - r2n * np.cos(nu2)
        ok &= np.abs(den) >= 1e-10 * np.maximum(r1n, r2n)
        e_t = (r2n - r1n) / den
        p_t = r1n * (1 + e_t * np.cos(nu1))

    out["penalty"][~ok] = 1e8
    g = ok                                    # geometry defined
    out["i_t"][g], out["raan_t"][g], out["nu1_t"][g], out["nu2_t"][g] = i_t[g], raan_t[g], nu1[g], nu2[g]
    out["arc_deg"][g] = np.rad2deg(np.mod(nu2[g] - nu1[g], TWO_PI))

    finite = g & np.isfinite(e_t) & np.isfinite(p_t)
    out["penalty"][g & ~finite] = 1e9
    with np.errstate(invalid="ignore"):
        shape_pen = (np.where(e_t < 0, 1e6 + 1e6 * np.abs(e_t), 0.0)
                     + np.where(e_t >= 1, 1e6 + 1e6 * (e_t - 1)**2, 0.0)
                     + np.where(p_t <= 0, 1e8 + 1e-2 * np.abs(p_t), 0.0))
    bad_shape = finite & (shape_pen > 0)
    out["penalty"][bad_shape] = shape_pen[bad_shape]
    out["e_t"][bad_shape] = e_t[bad_shape]

    ell = finite & ~bad_shape                 # proper ellipse
    a_t = np.where(ell, p_t / (1 - e_t**2), np.nan)
    rp_t = a_t * (1 - e_t)
    out["a_t"][ell], out["e_t"][ell], out["rp_t"][ell] = a_t[ell], e_t[ell], rp_t[ell]
    low = ell & (rp_t <= s.min_perihelion)
    out["penalty"][low] = 1e7 + (s.min_perihelion - rp_t[low])**2

    good = ell & ~low
    if np.any(good):
        _, v_t1 = kep2car(a_t[good], e_t[good], i_t[good], raan_t[good], X[good, 2], nu1[good], mu)
        _, v_t2 = kep2car(a_t[good], e_t[good], i_t[good], raan_t[good], X[good, 2], nu2[good], mu)
        out["v_t1"][good], out["v_t2"][good] = v_t1, v_t2
        dv1 = np.linalg.norm(v_t1 - v_earth[good], axis=1)
        dv2 = np.linalg.norm(v_target[good] - v_t2, axis=1)
        tof = time_of_flight(a_t[good], e_t[good], nu1[good], nu2[good], mu) / 86400
        out["dv1"][good], out["dv2"][good], out["dv"][good] = dv1, dv2, dv1 + dv2
        out["tof_days"][good] = tof
        J = dv1 + dv2 + s.weight_tof * tof
        out["J"][good] = J
        out["valid"][good] = np.isfinite(J) & (tof > 0)
    return out


def objective(X, cfg):
    """Cost for the optimizers: J if valid, otherwise 1e9 + penalty (like transferObjective.m)."""
    o = evaluate_transfers(X, cfg)
    J = np.where(o["valid"], o["J"], 1e9 + o["penalty"])
    return np.where(np.isfinite(J), J, 1e12)


def evaluate_transfer(x, cfg, method=""):
    """Single transfer as a namespace with scalar fields (vectors as arrays (3,))."""
    o = evaluate_transfers(np.asarray(x, float)[None, :], cfg)
    sol = NS(**{k: (v[0] if np.ndim(v) == 1 else v[0, :]) for k, v in o.items()})
    sol.valid = bool(sol.valid)
    sol.method = method
    return sol


def objective_and_gradient(x, cfg):
    """Cost and forward-difference gradient (step as fmincon) from a single vectorized call."""
    x = np.asarray(x, float)
    h = np.sqrt(np.finfo(float).eps) * np.maximum(1.0, np.abs(x))
    J = objective(np.vstack([x, x + np.diag(h)]), cfg)
    return J[0], (J[1:] - J[0]) / h


def run_scenario2(cfg):
    """Grid search, local optimization, differential evolution and multi-start (see runScenario2.m)."""
    s = cfg.sc2
    f = lambda x: float(objective(x, cfg)[0])                     # noqa: E731 scalar objective
    fg = lambda x: objective_and_gradient(x, cfg)                  # noqa: E731

    # 1) grid search, fully vectorized (loop order of the MATLAB code: argp_t fastest)
    vecs = [np.linspace(0, TWO_PI, n + 1)[:-1] for n in s.grid]
    N1, N2, N3 = np.meshgrid(*vecs, indexing="ij")
    X = np.stack([N1.ravel(), N2.ravel(), N3.ravel()], 1)
    rec = evaluate_transfers(X, cfg)
    J = np.where(rec["valid"], rec["J"], np.inf)
    r2 = NS(grid_vectors=vecs, grid_records=rec)
    r2.grid = evaluate_transfer(X[np.argmin(J)], cfg, "Grid search")

    # 2) local optimization from the best grid point (SLSQP ~ fmincon 'sqp', no bounds)
    history = []
    res = minimize(fg, r2.grid.x, method="SLSQP", jac=True, callback=lambda xk: history.append(f(xk)),
                   options=dict(maxiter=s.local.max_iterations, ftol=1e-14))
    r2.local = evaluate_transfer(res.x, cfg, "SLSQP")
    r2.local_history = np.array([r2.grid.J] + history)

    def polish(x0):
        return minimize(fg, x0, method="SLSQP", jac=True,
                        options=dict(maxiter=s.local.max_iterations, ftol=1e-14)).x

    # 3) differential evolution (counterpart of ga), several seeds, best one polished
    runs, best_de, best_hist = [], None, None
    for seed in s.de.seeds:
        hist = []
        res = differential_evolution(
            lambda Xt: objective(Xt.T, cfg), [(0, TWO_PI)] * 3, seed=seed,
            maxiter=s.de.generations, popsize=s.de.population // 3, tol=1e-8,
            polish=False, vectorized=True, updating="deferred",
            callback=lambda intermediate_result: hist.append(intermediate_result.fun))
        runs.append((seed, res.fun))
        if best_de is None or res.fun < best_de.fun:
            best_de, best_hist = res, hist
    r2.de_runs = np.array(runs)
    r2.de_history = np.array(best_hist)
    r2.de = evaluate_transfer(best_de.x, cfg, "Differential evolution")
    r2.de_refined = evaluate_transfer(polish(best_de.x), cfg, "DE + SLSQP")

    # 4) multi-start: SLSQP from the best grid point and from random points
    rng = np.random.default_rng(s.multistart.seed)
    starts = np.vstack([r2.grid.x, rng.uniform(-np.pi, 3 * np.pi, (s.multistart.points - 1, 3))])
    minima = [evaluate_transfer(polish(x0), cfg) for x0 in starts]
    valid = [m for m in minima if m.valid]
    best_ms = min(valid, key=lambda m: m.J)
    r2.multistart = evaluate_transfer(best_ms.x, cfg, "Multi-start")
    r2.local_minima = _distinct_minima(valid)

    candidates = [c for c in (r2.grid, r2.local, r2.de_refined, r2.multistart) if c.valid]
    r2.best = min(candidates, key=lambda c: c.J)
    r2.comparison = [r2.grid, r2.local, r2.de, r2.de_refined, r2.multistart]
    return r2


def _distinct_minima(sols):
    """Merge local minima that are the same transfer (angles shifted by 2*pi); count the starts."""
    table = {}
    for m in sols:
        key = (round(float(m.dv), 5), round(float(m.tof_days), 2))
        if key in table:
            table[key]["starts"] += 1
        else:
            table[key] = dict(dv=float(m.dv), tof_days=float(m.tof_days), arc_deg=float(m.arc_deg), starts=1)
    return sorted(table.values(), key=lambda r: r["dv"])
