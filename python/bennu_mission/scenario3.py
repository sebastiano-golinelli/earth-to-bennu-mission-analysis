"""Scenario 3: patched conics, Earth escape and capture at the target."""

from types import SimpleNamespace as NS

import numpy as np
from scipy.optimize import brentq

from .core import TWO_PI, ecliptic_to_equatorial, kep2car


def sphere_of_influence(distance_to_sun, mu_body, mu_sun):
    """Laplace sphere of influence d * (mu_body / mu_sun)^(2/5)."""
    return distance_to_sun * (mu_body / mu_sun) ** 0.4


def hyperbolic_leg(v_inf, rp, mu, r_body=0.0, r_soi=np.inf):
    """Planar hyperbola with excess speed v_inf and periapsis rp (a < 0)."""
    leg = NS(v_inf=v_inf, rp=rp, a=-mu / v_inf**2)
    leg.e = 1 - rp / leg.a
    leg.p = leg.a * (1 - leg.e**2)
    leg.v_esc = np.sqrt(2 * mu / rp)
    leg.v_circ = np.sqrt(mu / rp)
    leg.vp = np.sqrt(v_inf**2 + leg.v_esc**2)
    leg.impact_parameter = -leg.a * np.sqrt(leg.e**2 - 1)
    leg.turn_angle = 2 * np.arcsin(1 / leg.e)
    leg.nu_inf = np.arccos(-1 / leg.e)
    leg.feasible, leg.reason = True, ""
    if not (np.isfinite(v_inf) and v_inf > 0):
        leg.feasible, leg.reason = False, "invalid v_inf"
    elif rp <= r_body:
        leg.feasible, leg.reason = False, "periapsis below the surface"
    elif rp >= r_soi:
        leg.feasible, leg.reason = False, "periapsis outside the SOI"
    return leg


def hyperbolic_time_of_flight(a, e, mu, r):
    """Time [s] from periapsis to radius r on a hyperbola (a < 0)."""
    F = np.arccosh(max(1.0, (1 - r / a) / e))
    return (e * np.sinh(F) - F) / np.sqrt(mu / (-a) ** 3)


def capture_cost(v_inf, rp, mu):
    """Circular capture at the periapsis of the arrival hyperbola: (dv_circular, dv_marginal)."""
    vp_hyp = np.sqrt(v_inf**2 + 2 * mu / rp)
    return vp_hyp - np.sqrt(mu / rp), vp_hyp - np.sqrt(2 * mu / rp)


def capture_validity(rp, r_soi, mu_body, cfg, distance_to_sun):
    """r/rSOI, validity level, solar tide and SRP accelerations relative to the body's gravity."""
    ratio = rp / r_soi
    level = ("outside SOI" if ratio >= 1 else "low" if ratio > 0.5 else "medium" if ratio > 1 / 3 else "high")
    g_body = mu_body / rp**2
    tide = 2 * cfg.const.mu_sun * rp / distance_to_sun**3 / g_body
    s = cfg.sc3.srp
    a_srp = s.Cr * s.pressure_1au * (cfg.const.AU / distance_to_sun) ** 2 * s.area_to_mass / 1000
    return ratio, level, tide, a_srp / g_body


def escape_hyperbola_3d(v_inf, r_burn, v_park, mu, r_body=0.0):
    """Escape hyperbola through the burn point with outgoing asymptote parallel to v_inf.

    Same algorithm as escapeHyperbola3D.m; H.dv = |v_hyperbola - v_park| at the burn point.
    """
    v_inf, r_burn, v_park = (np.asarray(x, float) for x in (v_inf, r_burn, v_park))
    H = NS(a=np.nan, e=np.nan, i=np.nan, raan=np.nan, argp=np.nan, nu_burn=np.nan, nu_inf=np.nan,
           rp=np.nan, alpha=np.nan, r_burn=r_burn, v_burn=np.full(3, np.nan), dv=np.nan,
           feasible=False, reason="")
    v_mag, rb = np.linalg.norm(v_inf), np.linalg.norm(r_burn)
    if v_mag <= 0 or rb <= r_body:
        H.reason = "invalid input"
        return H
    u_inf, u_r = v_inf / v_mag, r_burn / rb
    a = -mu / v_mag**2
    alpha = np.arccos(np.clip(u_r @ u_inf, -1, 1))
    H.a, H.alpha = a, alpha

    e_max = 1 - rb / a
    if e_max <= 1.000001:
        H.reason = "no admissible eccentricity"
        return H

    def f(e):
        return rb - a * (1 - e**2) / (1 + e * np.cos(np.arccos(-1 / e) - alpha))

    e_grid = np.linspace(1 + 1e-6, e_max - 1e-9, 300)
    f_grid = f(e_grid)
    k = np.flatnonzero(f_grid[:-1] * f_grid[1:] < 0)
    if k.size == 0:
        H.reason = "incompatible geometry"
        return H
    e = brentq(f, e_grid[k[0]], e_grid[k[0] + 1], xtol=1e-15, rtol=4 * np.finfo(float).eps)
    H.e, H.nu_inf = e, np.arccos(-1 / e)
    H.nu_burn = H.nu_inf - alpha
    H.rp = a * (1 - e)

    h = np.cross(u_r, u_inf)
    if np.linalg.norm(h) < 1e-8:
        H.reason = "burn point aligned with v_inf"
        return H
    h /= np.linalg.norm(h)
    H.i = np.arccos(np.clip(h[2], -1, 1))
    n = np.cross([0.0, 0.0, 1.0], h)
    if np.linalg.norm(n) < 1e-12:
        n, H.raan = np.array([1.0, 0.0, 0.0]), 0.0
    else:
        n /= np.linalg.norm(n)
        H.raan = np.mod(np.arctan2(n[1], n[0]), TWO_PI)
    u = np.arctan2(np.cross(n, u_r) @ h, n @ u_r)      # argument of latitude of the burn point
    H.argp = np.mod(u - H.nu_burn, TWO_PI)

    r_check, H.v_burn = kep2car(a, e, H.i, H.raan, H.argp, H.nu_burn, mu)
    if np.linalg.norm(r_check - r_burn) > 1e-6 * rb:
        H.reason = "geometric inconsistency"
        return H
    H.dv = float(np.linalg.norm(H.v_burn - v_park))
    H.feasible = bool(H.rp > r_body)
    if not H.feasible:
        H.reason = "periapsis below the surface"
    return H


def run_scenario3(cfg, r1, r2):
    """Escape from the parking orbit of Scenario 1 towards the best transfer of Scenario 2, and capture."""
    c, s, tgt = cfg.const, cfg.sc3, cfg.sc2.target
    park = r1.final
    r3 = NS()
    r3.parking = NS(**vars(park), rp=park.a * (1 - park.e), ra=park.a * (1 + park.e),
                    vp=np.sqrt(c.mu_earth * (2 / (park.a * (1 - park.e)) - 1 / park.a)))
    best = r2.best
    r3.v_inf_dep_ecliptic = best.v_t1 - best.v_earth
    r3.v_inf_arr = best.v_t2 - best.v_target
    r3.v_inf_dep = float(np.linalg.norm(r3.v_inf_dep_ecliptic))
    r3.v_inf_arr_mag = float(np.linalg.norm(r3.v_inf_arr))

    # spheres of influence
    r3.soi_earth = sphere_of_influence(cfg.sc2.earth.a, c.mu_earth, c.mu_sun)
    r3.soi_target = sphere_of_influence(tgt.orbit.a, tgt.mu, c.mu_sun)
    d1, d2 = np.linalg.norm(best.r1), np.linalg.norm(best.r2)
    if s.soi_mode == "patch":
        soi_e, soi_t = sphere_of_influence(d1, c.mu_earth, c.mu_sun), sphere_of_influence(d2, tgt.mu, c.mu_sun)
    else:
        soi_e, soi_t = r3.soi_earth, r3.soi_target

    # coplanar escape: tangential burn at parking periapsis
    dep = hyperbolic_leg(r3.v_inf_dep, r3.parking.rp, c.mu_earth, c.r_earth, soi_e)
    if not dep.feasible:
        raise RuntimeError(f"escape hyperbola not feasible: {dep.reason}")
    dep.dv = dep.vp - r3.parking.vp
    dep.time_in_soi = hyperbolic_time_of_flight(dep.a, dep.e, c.mu_earth, soi_e)
    r3.escape = dep

    # capture at the target
    h = s.capture_altitudes
    rows = []
    for hk in h:
        r_cap = tgt.radius + hk
        leg = hyperbolic_leg(r3.v_inf_arr_mag, r_cap, tgt.mu, tgt.radius, soi_t)
        dv_circ, dv_marg = capture_cost(r3.v_inf_arr_mag, r_cap, tgt.mu)
        ratio, level, tide, srp = capture_validity(r_cap, soi_t, tgt.mu, cfg, d2)
        rows.append(dict(altitude_km=hk, radius_km=r_cap, hyperbola_ecc=leg.e, dv_circular=dv_circ,
                         dv_marginal=dv_marg, radius_over_soi=ratio, validity=level,
                         tide_ratio=tide, srp_ratio=srp))
    r3.capture_table = rows
    k_nom = int(np.argmin(np.abs(h - s.nominal_capture_altitude)))
    r3.capture = NS(**rows[k_nom])
    leg = hyperbolic_leg(r3.v_inf_arr_mag, tgt.radius + h[k_nom], tgt.mu)
    r3.capture.time_in_soi = hyperbolic_time_of_flight(leg.a, leg.e, tgt.mu, soi_t)
    r3.capture.dv = r3.capture.dv_circular

    # non-coplanar escape with the burn point optimized along the parking orbit
    r3.v_inf_dep_equatorial = ecliptic_to_equatorial(r3.v_inf_dep_ecliptic, c.obliquity)
    nu = np.linspace(0, TWO_PI, s.burn_sweep_points)
    rb, vb = kep2car(park.a, park.e, park.i, park.raan, park.argp, nu, c.mu_earth)
    dv = np.full(nu.size, np.nan)
    for k in range(nu.size):
        Hk = escape_hyperbola_3d(r3.v_inf_dep_equatorial, rb[k], vb[k], c.mu_earth, c.r_earth)
        if Hk.feasible:
            dv[k] = Hk.dv
    if np.all(np.isnan(dv)):
        raise RuntimeError("no feasible non-coplanar escape along the parking orbit")
    k_opt = int(np.nanargmin(dv))
    r3.escape_3d = escape_hyperbola_3d(r3.v_inf_dep_equatorial, rb[k_opt], vb[k_opt], c.mu_earth, c.r_earth)
    r3.escape_3d.nu_burn_on_parking = nu[k_opt]
    r3.burn_sweep = NS(nu=nu, dv=dv)

    r3.dv_coplanar = r3.escape.dv + r3.capture.dv
    r3.dv_noncoplanar = r3.escape_3d.dv + r3.capture.dv
    return r3
