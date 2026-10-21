"""Computations behind the web app, independent of Streamlit so that they can be tested directly.

Everything is built on the bennu_mission package, used unchanged.
"""

import sys
from pathlib import Path
from types import SimpleNamespace as NS
from typing import NamedTuple

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from bennu_mission import mission_config  # noqa: E402
from bennu_mission.core import TWO_PI, Orbit, ecliptic_to_equatorial  # noqa: E402
from bennu_mission.scenario1 import run_scenario1  # noqa: E402
from bennu_mission.scenario2 import evaluate_transfer, evaluate_transfers, objective_and_gradient  # noqa: E402
from bennu_mission.scenario3 import run_scenario3  # noqa: E402

D2R = np.pi / 180
DEFAULT = mission_config()
AU = DEFAULT.const.AU
BENNU = DEFAULT.sc2.target.orbit


class Params(NamedTuple):
    """Everything the user can change; hashable, so results can be cached."""

    a_au: float = BENNU.a / AU
    e: float = BENNU.e
    i_deg: float = BENNU.i / D2R
    raan_deg: float = BENNU.raan / D2R
    argp_deg: float = BENNU.argp / D2R
    min_perihelion_au: float = 0.3
    allow_long_arc: bool = True
    tof_penalty: float = 0.0      # [m/s per day of flight]


def is_bennu(p: Params):
    ref = Params()
    return bool(np.allclose(p[:5], ref[:5], rtol=1e-9, atol=1e-9))


def make_cfg(p: Params):
    cfg = mission_config()
    cfg.verbose = False
    cfg.plot.enable = False
    s = cfg.sc2
    s.target.orbit = Orbit(a=p.a_au * AU, e=p.e, i=p.i_deg * D2R, raan=p.raan_deg * D2R, argp=p.argp_deg * D2R)
    if not is_bennu(p):
        s.target.name = "Custom target"
    s.min_perihelion = p.min_perihelion_au * AU
    s.allow_long_arc = p.allow_long_arc
    s.weight_tof = p.tof_penalty / 1000
    return cfg


def refine_candidates(X, J, cfg, count=12, separation=20 * D2R):
    """SLSQP from the best grid points that are at least `separation` apart in every angle.

    Refining only the best grid point can stop in a local minimum when the target changes; a dozen
    distinct starting points recover the optimum of the package's full global search (differential
    evolution + 200-point multi-start) at a fraction of the cost.
    """
    picks = []
    for n in np.argsort(J):
        if not np.isfinite(J[n]) or len(picks) == count:
            break
        if not picks or np.all(np.abs((X[picks] - X[n] + np.pi) % TWO_PI - np.pi).max(axis=1) > separation):
            picks.append(n)
    sols = [evaluate_transfer(X[picks[0]], cfg, "Grid search")]
    for n in picks:
        x = minimize(lambda x: objective_and_gradient(x, cfg), X[n], method="SLSQP", jac=True,
                     options=dict(maxiter=cfg.sc2.local.max_iterations, ftol=1e-14)).x
        sols.append(evaluate_transfer(x, cfg, "Grid search + SLSQP"))
    return min((s for s in sols if s.valid), key=lambda s: s.J)


def optimize_transfer(p: Params):
    """Vectorized grid search over departure point, arrival point and transfer orientation, then refinement.

    Returns None when no grid point gives a valid transfer.
    """
    cfg = make_cfg(p)
    vecs = [np.linspace(0, TWO_PI, n + 1)[:-1] for n in cfg.sc2.grid]
    N1, N2, N3 = np.meshgrid(*vecs, indexing="ij")
    X = np.stack([N1.ravel(), N2.ravel(), N3.ravel()], 1)
    rec = evaluate_transfers(X, cfg)
    J = np.where(rec["valid"], rec["J"], np.inf)
    if not np.isfinite(J).any():
        return None
    best = refine_candidates(X, J, cfg)
    cost = J.reshape(cfg.sc2.grid).min(axis=2)              # best orientation for each departure/arrival pair
    cost[np.isinf(cost)] = np.nan
    v = rec["valid"]
    return dict(best=best, cost=cost, vecs=[np.rad2deg(x) for x in vecs], n_valid=int(v.sum()), n_total=int(v.size),
                tof=rec["tof_days"][v], dv=rec["dv"][v], arc=rec["arc_deg"][v])


def close_mission(p: Params, best):
    """Scenarios 1 and 3 around the heliocentric transfer `best` (same chain and budget as run_mission)."""
    cfg = make_cfg(p)
    v_inf = ecliptic_to_equatorial(best.v_t1 - best.v_earth, cfg.const.obliquity)
    r1 = run_scenario1(cfg, v_inf)
    r3 = run_scenario3(cfg, r1, NS(best=best))
    days = r1.best.time / 86400 + r3.escape.time_in_soi / 86400 + best.tof_days + r3.capture.time_in_soi / 86400
    return dict(
        launch_raan=float(np.rad2deg(r1.initial.raan)), sc1_name=r1.best.name, sc1_dv=float(r1.best.dv_total),
        sc1_hours=float(r1.best.time / 3600), sweep_raan=np.asarray(r1.launch_sweep.raan_deg),
        sweep_dv=np.asarray(r1.launch_sweep.dv),
        sc1_table=[{"Strategy": c.name, "Δv [km/s]": round(float(c.dv_total), 4),
                    "Time [h]": round(float(np.min(c.dt_total) / 3600), 2)}
                   for c in sorted(r1.cases, key=lambda c: c.dv_total)],
        escape_2d=float(r3.escape.dv), escape_3d=float(r3.escape_3d.dv), v_inf_dep=float(r3.v_inf_dep),
        capture=float(r3.capture.dv), capture_alt=float(r3.capture.altitude_km), v_inf_arr=float(r3.v_inf_arr_mag),
        soi_target=float(r3.soi_target), escape_days=float(r3.escape.time_in_soi / 86400),
        total=float(r1.best.dv_total + r3.dv_noncoplanar), days=float(days), tof=float(best.tof_days))
