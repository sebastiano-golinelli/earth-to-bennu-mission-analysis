"""Impulsive maneuvers for the geocentric transfer (scalar inputs)."""

from dataclasses import dataclass, field

import numpy as np

from .core import TWO_PI, Orbit, time_of_flight


def _clip(x):
    return float(np.clip(x, -1.0, 1.0))


def bitangent_transfer(a_i, e_i, a_f, e_f, kind, argp_i, mu):
    """Two-impulse bitangent transfer between coaxial ellipses.

    kind = 'pa', 'ap', 'pp' or 'aa': departure and arrival apses. For 'pp'
    and 'aa' the apse line of the final orbit is rotated by pi.
    Returns (dv1, dv2, dt, argp_f, a_t, e_t); dv1 and dv2 are signed.
    """
    if kind == "pa":
        r1, r2, d_argp = a_i * (1 - e_i), a_f * (1 + e_f), 0.0
    elif kind == "ap":
        r1, r2, d_argp = a_i * (1 + e_i), a_f * (1 - e_f), 0.0
    elif kind == "pp":
        r1, r2, d_argp = a_i * (1 - e_i), a_f * (1 - e_f), np.pi
    elif kind == "aa":
        r1, r2, d_argp = a_i * (1 + e_i), a_f * (1 + e_f), np.pi
    else:
        raise ValueError("bitangent_transfer: kind must be pa, ap, pp or aa")
    a_t = (r1 + r2) / 2
    e_t = abs(r2 - r1) / (r1 + r2)
    dv1 = np.sqrt(mu) * (np.sqrt(2 / r1 - 1 / a_t) - np.sqrt(2 / r1 - 1 / a_i))
    dv2 = np.sqrt(mu) * (np.sqrt(2 / r2 - 1 / a_f) - np.sqrt(2 / r2 - 1 / a_t))
    argp_f = np.mod(argp_i + d_argp, TWO_PI)
    dt = np.pi * np.sqrt(a_t**3 / mu)
    return dv1, dv2, dt, argp_f, a_t, e_t


def plane_change(a, e, i_i, raan_i, argp_i, i_f, raan_f, mu):
    """Single-impulse plane change at the node with the lower speed.

    Returns (dv, argp_f, nu): impulse, new argument of periapsis and true
    anomaly of the maneuver.
    """
    d_raan = np.mod(raan_f - raan_i + np.pi, TWO_PI) - np.pi
    d_i = i_f - i_i
    alpha = np.arccos(_clip(np.cos(i_i) * np.cos(i_f) + np.sin(i_i) * np.sin(i_f) * np.cos(d_raan)))
    if alpha < 1e-12:
        return 0.0, argp_i, 0.0
    sa = np.sin(alpha)
    if abs(sa) < 1e-12 or abs(np.sin(i_i)) < 1e-12 or abs(np.sin(i_f)) < 1e-12:
        raise ValueError("plane_change: degenerate geometry (alpha = 0/pi or equatorial orbit)")

    # argument of latitude of the node on the initial (u_i) and final (u_f) plane
    s_ui = _clip(np.sin(abs(d_raan)) * np.sin(i_f) / sa)
    s_uf = _clip(np.sin(abs(d_raan)) * np.sin(i_i) / sa)
    if d_i >= 0:
        c_ui = _clip((-np.cos(i_f) + np.cos(alpha) * np.cos(i_i)) / (sa * np.sin(i_i)))
        c_uf = _clip((np.cos(i_i) - np.cos(alpha) * np.cos(i_f)) / (sa * np.sin(i_f)))
    else:
        c_ui = _clip((np.cos(i_f) - np.cos(alpha) * np.cos(i_i)) / (sa * np.sin(i_i)))
        c_uf = _clip((-np.cos(i_i) + np.cos(alpha) * np.cos(i_f)) / (sa * np.sin(i_f)))
    u_i = np.mod(np.arctan2(s_ui, c_ui), TWO_PI)
    u_f = np.mod(np.arctan2(s_uf, c_uf), TWO_PI)

    if (d_raan >= 0) == (d_i >= 0):
        nu = np.mod(u_i - argp_i, TWO_PI)
        argp_f = np.mod(u_f - nu, TWO_PI)
    else:
        nu = np.mod(TWO_PI - u_i - argp_i, TWO_PI)
        argp_f = np.mod(TWO_PI - u_f - nu, TWO_PI)

    if np.cos(nu) > 0:                      # use the node on the apoapsis side
        nu = np.mod(nu + np.pi, TWO_PI)

    p = a * (1 - e**2)
    v_theta = np.sqrt(mu / p) * (1 + e * np.cos(nu))
    return 2 * v_theta * np.sin(alpha / 2), argp_f, nu


def periapsis_change(a, e, argp_i, argp_f, mu):
    """Rotation of the apse line in the orbit plane.

    Returns (dv, nu_i, nu_f): nu_i and nu_f (arrays of 2) are the true
    anomalies of the two possible maneuver points on the initial and final orbit.
    """
    if e < 1e-10:
        raise ValueError("periapsis_change: argument of periapsis undefined for a circular orbit")
    if e >= 1:
        raise ValueError("periapsis_change: only elliptic orbits (e < 1) are supported")
    d = np.mod(argp_f - argp_i, TWO_PI)
    if d < 1e-12 or abs(d - TWO_PI) < 1e-12:
        return 0.0, np.array([0.0, np.pi]), np.array([0.0, np.pi])
    nu_i = np.mod(np.array([d / 2, np.pi + d / 2]), TWO_PI)
    nu_f = np.mod(np.array([TWO_PI - d / 2, np.pi - d / 2]), TWO_PI)
    p = a * (1 - e**2)
    return 2 * e * np.sqrt(mu / p) * abs(np.sin(d / 2)), nu_i, nu_f


@dataclass
class BitangentPlaneChange:
    """Result of bitangent_with_plane_change (angles in rad, times in s)."""

    dv1: float = np.nan
    dv2: float = np.nan
    dv_plane: float = np.nan
    argp_before: float = np.nan
    argp_after: float = np.nan
    argp_final: float = np.nan
    a_t: float = np.nan
    e_t: float = np.nan
    nu_node: float = np.nan
    dt_coast: np.ndarray = field(default_factory=lambda: np.full(2, np.nan))
    feasible: bool = False


def _on_arc(nu, nu0, nu1):
    nu, nu0, nu1 = np.mod(nu, TWO_PI), np.mod(nu0, TWO_PI), np.mod(nu1, TWO_PI)
    if nu0 <= nu1:
        return nu0 <= nu <= nu1
    return nu >= nu0 or nu <= nu1


def bitangent_with_plane_change(o_i: Orbit, o_f: Orbit, kind, mu):
    """Bitangent transfer with the plane change moved onto the transfer arc.

    The plane change is done at the node between the transfer plane and the
    final plane. Assumes that for pa/pp/aa the departure apsis is lower than
    the arrival apsis, and higher for ap.
    """
    if kind == "pa":
        r_dep, r_arr = o_i.a * (1 - o_i.e), o_f.a * (1 + o_f.e)
        nu_dep, nu_arr, argp_before, d_final = 0.0, np.pi, o_i.argp, 0.0
    elif kind == "ap":
        r_dep, r_arr = o_i.a * (1 + o_i.e), o_f.a * (1 - o_f.e)
        nu_dep, nu_arr, argp_before, d_final = np.pi, TWO_PI, o_i.argp, 0.0
    elif kind == "pp":
        r_dep, r_arr = o_i.a * (1 - o_i.e), o_f.a * (1 - o_f.e)
        nu_dep, nu_arr, argp_before, d_final = 0.0, np.pi, o_i.argp, np.pi
    elif kind == "aa":
        r_dep, r_arr = o_i.a * (1 + o_i.e), o_f.a * (1 + o_f.e)
        nu_dep, nu_arr, argp_before, d_final = 0.0, np.pi, np.mod(o_i.argp + np.pi, TWO_PI), 0.0
    else:
        raise ValueError("bitangent_with_plane_change: kind must be pa, ap, pp or aa")
    if (kind == "ap") != (r_dep > r_arr):
        raise ValueError(f"bitangent_with_plane_change: unsupported apsis ordering for type {kind}")

    m = BitangentPlaneChange(argp_before=argp_before)
    m.a_t = (r_dep + r_arr) / 2
    m.e_t = abs(r_arr - r_dep) / (r_arr + r_dep)
    m.dv1 = np.sqrt(mu) * (np.sqrt(2 / r_dep - 1 / m.a_t) - np.sqrt(2 / r_dep - 1 / o_i.a))
    m.dv2 = np.sqrt(mu) * (np.sqrt(2 / r_arr - 1 / o_f.a) - np.sqrt(2 / r_arr - 1 / m.a_t))
    m.dv_plane, m.argp_after, m.nu_node = plane_change(m.a_t, m.e_t, o_i.i, o_i.raan, argp_before,
                                                       o_f.i, o_f.raan, mu)
    if not _on_arc(m.nu_node, nu_dep, nu_arr):
        nu_other = np.mod(m.nu_node + np.pi, TWO_PI)       # try the other node
        if not _on_arc(nu_other, nu_dep, nu_arr):
            return m
        m.nu_node = nu_other
        p = m.a_t * (1 - m.e_t**2)
        v_node = np.sqrt(mu / p) * (1 + m.e_t * np.cos(m.nu_node))
        alpha = np.arccos(_clip(np.cos(o_i.i) * np.cos(o_f.i)
                                + np.sin(o_i.i) * np.sin(o_f.i) * np.cos(o_f.raan - o_i.raan)))
        m.dv_plane = 2 * v_node * np.sin(alpha / 2)

    m.feasible = True
    m.dt_coast = np.array([time_of_flight(m.a_t, m.e_t, nu_dep, m.nu_node, mu),
                           time_of_flight(m.a_t, m.e_t, m.nu_node, nu_arr, mu)])
    m.argp_final = np.mod(m.argp_after + d_final, TWO_PI)
    return m
