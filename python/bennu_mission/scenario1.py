"""Scenario 1: geocentric transfer from the launch orbit to the parking orbit."""

from dataclasses import dataclass, field, replace
from types import SimpleNamespace as NS

import numpy as np

from .core import TWO_PI, Orbit, apsis_anomaly, car2kep, kep2car, time_of_flight
from .maneuvers import bitangent_transfer, bitangent_with_plane_change, periapsis_change, plane_change


@dataclass
class Strategy:
    """Cost and duration of one strategy.

    dv = [dv1 bitangent, dv2 bitangent, dv plane, dv periapsis] [km/s];
    dt_total = total time [s] for the two possible periapsis-change points.
    """

    name: str
    strategy: str
    kind: str
    feasible: bool = True
    dv: np.ndarray = field(default_factory=lambda: np.full(4, np.nan))
    dv_total: float = np.nan
    dt_total: np.ndarray = field(default_factory=lambda: np.full(2, np.nan))
    geo: object = None

    @property
    def time(self):
        return float(np.min(self.dt_total))


def evaluate_strategy(strategy, kind, o_i: Orbit, o_f: Orbit, mu):
    """Transfer from o_i (nu = current position) to o_f (nu = arrival point).

    strategy: 'SPW' shape -> plane -> periapsis, 'PSW', 'PWS', or 'S(P)W'
    (shape with the plane change inside the transfer arc -> periapsis);
    kind: bitangent type 'ap', 'pa', 'pp', 'aa'.
    """
    c = Strategy(name=f"{strategy} {kind}", strategy=strategy, kind=kind)
    nu_dep = apsis_anomaly(kind[0])
    nu_arr = apsis_anomaly(kind[1])

    def tof(a, e, nu1, nu2):
        return time_of_flight(a, e, nu1, nu2, mu)

    if strategy == "SPW":
        dt1 = tof(o_i.a, o_i.e, o_i.nu, nu_dep)
        dv1, dv2, dt2, argp_s, _, _ = bitangent_transfer(o_i.a, o_i.e, o_f.a, o_f.e, kind, o_i.argp, mu)
        dvp, argp_p, nu_p = plane_change(o_f.a, o_f.e, o_i.i, o_i.raan, argp_s, o_f.i, o_f.raan, mu)
        dt3 = tof(o_f.a, o_f.e, nu_arr, nu_p)
        dvw, nu_wi, nu_wf = periapsis_change(o_f.a, o_f.e, argp_p, o_f.argp, mu)
        dt = dt1 + dt2 + dt3 + tof(o_f.a, o_f.e, nu_p, nu_wi) + tof(o_f.a, o_f.e, nu_wf, o_f.nu)
        c.geo = NS(argp_shape=argp_s, argp_plane=argp_p)

    elif strategy == "PSW":
        dvp, argp_p, nu_p = plane_change(o_i.a, o_i.e, o_i.i, o_i.raan, o_i.argp, o_f.i, o_f.raan, mu)
        dt1 = tof(o_i.a, o_i.e, o_i.nu, nu_p)
        dt2 = tof(o_i.a, o_i.e, nu_p, nu_dep)
        dv1, dv2, dt3, argp_s, _, _ = bitangent_transfer(o_i.a, o_i.e, o_f.a, o_f.e, kind, argp_p, mu)
        dvw, nu_wi, nu_wf = periapsis_change(o_f.a, o_f.e, argp_s, o_f.argp, mu)
        dt = dt1 + dt2 + dt3 + tof(o_f.a, o_f.e, nu_arr, nu_wi) + tof(o_f.a, o_f.e, nu_wf, o_f.nu)
        c.geo = NS(argp_plane=argp_p, argp_shape=argp_s)

    elif strategy == "PWS":
        dvp, argp_p, nu_p = plane_change(o_i.a, o_i.e, o_i.i, o_i.raan, o_i.argp, o_f.i, o_f.raan, mu)
        # periapsis target before the bitangent (pp and aa rotate the apse line by pi)
        argp_w = np.mod(o_f.argp - np.pi, TWO_PI) if kind[0] == kind[1] else o_f.argp
        dt1 = tof(o_i.a, o_i.e, o_i.nu, nu_p)
        dvw, nu_wi, nu_wf = periapsis_change(o_i.a, o_i.e, argp_p, argp_w, mu)
        dt2 = tof(o_i.a, o_i.e, nu_p, nu_wi)
        dv1, dv2, dt3, _, _, _ = bitangent_transfer(o_i.a, o_i.e, o_f.a, o_f.e, kind, argp_w, mu)
        dt4 = tof(o_i.a, o_i.e, nu_wf, nu_dep)
        dt5 = tof(o_f.a, o_f.e, nu_arr, o_f.nu)
        dt = dt1 + dt2 + dt4 + dt3 + dt5
        c.geo = NS(argp_plane=argp_p, argp_periapsis=argp_w)

    elif strategy == "S(P)W":
        dt1 = tof(o_i.a, o_i.e, o_i.nu, nu_dep)
        m = bitangent_with_plane_change(o_i, o_f, kind, mu)
        if not m.feasible:
            c.feasible = False
            return c
        dv1, dv2, dvp = m.dv1, m.dv2, m.dv_plane
        dvw, nu_wi, nu_wf = periapsis_change(o_f.a, o_f.e, m.argp_final, o_f.argp, mu)
        dt = dt1 + m.dt_coast.sum() + tof(o_f.a, o_f.e, nu_arr, nu_wi) + tof(o_f.a, o_f.e, nu_wf, o_f.nu)
        c.geo = m

    else:
        raise ValueError(f"evaluate_strategy: unknown strategy {strategy}")

    c.dv = np.array([dv1, dv2, dvp, dvw], dtype=float)
    c.dv_total = float(np.sum(np.abs(c.dv)))
    c.dt_total = np.broadcast_to(np.asarray(dt, float), (2,)).copy()
    return c


def design_parking_orbit(v_inf, rp, ra, ref: Orbit, mu):
    """Parking orbit whose tangential periapsis burn escapes along v_inf (equatorial frame).

    The plane contains v_inf and has the inclination of ref if possible
    (otherwise the declination of v_inf); of the two such planes, the one
    closer to the plane of ref is chosen. The periapsis precedes the
    asymptote by the true anomaly of the asymptote. Arrival point: nu = 0.
    """
    u = np.asarray(v_inf, float) / np.linalg.norm(v_inf)
    dec = np.arcsin(u[2])
    ra0 = np.arctan2(u[1], u[0])
    i = max(ref.i, abs(dec))
    k = float(np.clip(np.tan(dec) / np.tan(i), -1, 1))
    raan_cand = np.mod(np.array([ra0 - np.arcsin(k), ra0 - np.pi + np.arcsin(k)]), TWO_PI)
    plane_angle = np.arccos(np.clip(np.cos(ref.i) * np.cos(i)
                                    + np.sin(ref.i) * np.sin(i) * np.cos(raan_cand - ref.raan), -1, 1))
    raan = raan_cand[np.argmin(plane_angle)]

    n = np.array([np.cos(raan), np.sin(raan), 0.0])
    h = np.array([np.sin(raan) * np.sin(i), -np.cos(raan) * np.sin(i), np.cos(i)])
    u_inf = np.arctan2(u @ np.cross(h, n), u @ n)
    e_hyp = 1 + rp * np.linalg.norm(v_inf)**2 / mu
    nu_inf = np.arccos(-1 / e_hyp)
    park = Orbit(a=(rp + ra) / 2, e=(ra - rp) / (ra + rp), i=i, raan=raan,
                 argp=np.mod(u_inf - nu_inf, TWO_PI), nu=0.0)
    return park, NS(right_ascension=ra0, declination=dec)


def _parking_orbit(cfg, v_inf, o_i):
    p = cfg.sc1.parking
    if p.mode == "designed":
        park, _ = design_parking_orbit(v_inf, cfg.const.r_earth + p.periapsis_altitude,
                                       cfg.const.r_earth + p.apoapsis_altitude, o_i, cfg.const.mu_earth)
        return park
    if p.mode == "given":
        return car2kep(p.state_r, p.state_v, cfg.const.mu_earth)
    raise ValueError("parking.mode must be 'designed' or 'given'")


def _evaluate_all(cfg, o_i, o_f):
    return [evaluate_strategy(s, t, o_i, o_f, cfg.const.mu_earth)
            for s in cfg.sc1.strategies for t in cfg.sc1.types]


def _select(cases, selection):
    if selection == "minDeltaV":
        return int(np.nanargmin([c.dv_total for c in cases]))
    if selection == "minTime":
        return int(np.nanargmin([np.min(c.dt_total) for c in cases]))
    names = [c.name for c in cases]
    if selection not in names:
        raise ValueError(f"unknown selection {selection}")
    return names.index(selection)


def run_scenario1(cfg, v_inf):
    """All strategies from the launch orbit to the parking orbit (see runScenario1.m)."""
    o_i = cfg.sc1.initial
    r1 = NS(launch_sweep=None)
    if cfg.sc1.launch_raan == "optimize":
        raan_deg = np.arange(0, 360, cfg.sc1.launch_raan_step, dtype=float)
        dv = np.empty(raan_deg.size)
        time = np.empty(raan_deg.size)
        name = []
        for k, raan in enumerate(raan_deg):
            o_k = replace(o_i, raan=np.deg2rad(raan))
            cases = _evaluate_all(cfg, o_k, _parking_orbit(cfg, v_inf, o_k))
            sel = cases[_select(cases, cfg.sc1.selection)]
            dv[k], time[k] = sel.dv_total, sel.time / 3600
            name.append(sel.name)
        r1.launch_sweep = NS(raan_deg=raan_deg, dv=dv, time_h=time, strategy=name)
        o_i = replace(o_i, raan=np.deg2rad(raan_deg[np.argmin(dv)]))
    else:
        o_i = replace(o_i, raan=np.deg2rad(cfg.sc1.launch_raan))

    r1.initial = o_i
    r1.initial_state = kep2car(o_i.a, o_i.e, o_i.i, o_i.raan, o_i.argp, o_i.nu, cfg.const.mu_earth)
    r1.final = _parking_orbit(cfg, v_inf, o_i)
    if cfg.sc1.parking.mode == "designed":
        _, r1.asymptote = design_parking_orbit(v_inf, cfg.const.r_earth + cfg.sc1.parking.periapsis_altitude,
                                               cfg.const.r_earth + cfg.sc1.parking.apoapsis_altitude,
                                               o_i, cfg.const.mu_earth)
    r1.cases = _evaluate_all(cfg, o_i, r1.final)
    r1.selected = _select(r1.cases, cfg.sc1.selection)
    r1.best = r1.cases[r1.selected]
    return r1
