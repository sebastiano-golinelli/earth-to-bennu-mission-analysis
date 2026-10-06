"""Basic orbital mechanics: Keplerian elements, frames, time of flight.

Units: km, km/s, s, rad. Vectors are numpy arrays with the three components
on the last axis, so every function also works on stacks of states
(shape (N, 3)) without loops.
"""

from dataclasses import dataclass

import numpy as np

TWO_PI = 2.0 * np.pi


@dataclass
class Orbit:
    """Keplerian elements; nu is the true anomaly of the point of interest."""

    a: float
    e: float
    i: float
    raan: float
    argp: float
    nu: float = 0.0


def _perifocal_axes(raan, i, argp):
    """Inertial unit vectors of the perifocal axes p (to periapsis), q and w (normal).

    They are the columns of R3(raan) R1(i) R3(argp), written in closed form
    to avoid building small matrices (much faster for single points).
    """
    cO, sO = np.cos(raan), np.sin(raan)
    ci, si = np.cos(i), np.sin(i)
    cw, sw = np.cos(argp), np.sin(argp)
    p = np.stack(np.broadcast_arrays(cO * cw - sO * ci * sw, sO * cw + cO * ci * sw, si * sw), -1)
    q = np.stack(np.broadcast_arrays(-cO * sw - sO * ci * cw, -sO * sw + cO * ci * cw, si * cw), -1)
    w = np.stack(np.broadcast_arrays(sO * si, -cO * si, ci), -1)
    return p, q, w


def perifocal_to_inertial(raan, i, argp):
    """Rotation matrix R3(raan) R1(i) R3(argp), shape (..., 3, 3)."""
    return np.stack(_perifocal_axes(raan, i, argp), -1)


def kep2car(a, e, i, raan, argp, nu, mu):
    """Keplerian elements -> inertial position [km] and velocity [km/s].

    Works for ellipses and hyperbolas (a < 0). All inputs broadcast; the
    outputs have shape (..., 3).
    """
    a, e, nu = np.asarray(a, float), np.asarray(e, float), np.asarray(nu, float)
    p_axis, q_axis, _ = _perifocal_axes(raan, i, argp)
    p = a * (1.0 - e**2)
    c, s = np.cos(nu)[..., None], np.sin(nu)[..., None]
    r = (p / (1.0 + e * np.cos(nu)))[..., None] * (c * p_axis + s * q_axis)
    v = np.sqrt(mu / p)[..., None] * (-s * p_axis + (e[..., None] + c) * q_axis)
    return r, v


def car2kep(r, v, mu):
    """Inertial position and velocity -> Orbit.

    Raises ValueError for equatorial or circular orbits, where the RAAN or
    the argument of periapsis is undefined.
    """
    r, v = np.asarray(r, float), np.asarray(v, float)
    rn = np.linalg.norm(r)
    h = np.cross(r, v)
    n = np.cross([0.0, 0.0, 1.0], h)
    e_vec = np.cross(v, h) / mu - r / rn
    e = np.linalg.norm(e_vec)
    if np.linalg.norm(n) < 1e-12 * np.linalg.norm(h) or e < 1e-12:
        raise ValueError("car2kep: equatorial or circular orbit, RAAN or argument of periapsis undefined")
    n = n / np.linalg.norm(n)

    a = 1.0 / (2.0 / rn - v @ v / mu)
    i = np.arccos(np.clip(h[2] / np.linalg.norm(h), -1, 1))
    raan = np.arccos(np.clip(n[0], -1, 1))
    if n[1] < 0:
        raan = TWO_PI - raan
    argp = np.arccos(np.clip(n @ e_vec / e, -1, 1))
    if e_vec[2] < 0:
        argp = TWO_PI - argp
    nu = np.arccos(np.clip(r @ e_vec / (rn * e), -1, 1))
    if r @ v < 0:
        nu = TWO_PI - nu
    return Orbit(float(a), float(e), float(i), float(raan), float(argp), float(nu))


def time_of_flight(a, e, nu1, nu2, mu):
    """Time [s] from true anomaly nu1 to nu2 on an ellipse, in the direction of motion.

    Returns 0 <= dt < period (0 when nu1 == nu2). Inputs broadcast.
    """
    e = np.asarray(e, float)
    if np.any(e >= 1):
        raise ValueError("time_of_flight: only elliptic orbits (e < 1) are supported")
    nu1 = np.mod(nu1, TWO_PI)
    nu2 = np.mod(nu2, TWO_PI)
    E1 = np.mod(2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu1 / 2), np.sqrt(1 + e) * np.cos(nu1 / 2)), TWO_PI)
    E2 = np.mod(2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu2 / 2), np.sqrt(1 + e) * np.cos(nu2 / 2)), TWO_PI)
    M1 = E1 - e * np.sin(E1)
    M2 = E2 - e * np.sin(E2)
    return np.mod(M2 - M1, TWO_PI) / np.sqrt(mu / np.asarray(a, float)**3)


def conic_points(a, e, i, raan, argp, nu):
    """Positions (N, 3) along a conic for an array of true anomalies (for plots)."""
    r, _ = kep2car(a, e, i, raan, argp, np.asarray(nu, float), 1.0)
    return r


def apsis_anomaly(letter):
    """'p' (periapsis) -> 0, 'a' (apoapsis) -> pi."""
    if letter == "p":
        return 0.0
    if letter == "a":
        return np.pi
    raise ValueError("apsis_anomaly: use 'p' or 'a'")


def ecliptic_to_equatorial(v, obliquity):
    """Rotate a vector from the ecliptic to the equatorial frame (about the x axis)."""
    c, s = np.cos(obliquity), np.sin(obliquity)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]]) @ np.asarray(v, float)
