"""Element conversions, time of flight and maneuver geometry (same checks as the MATLAB tests)."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from bennu_mission.core import car2kep, kep2car, time_of_flight
from bennu_mission.maneuvers import bitangent_transfer, periapsis_change, plane_change

MU = 398600.4418


def angle_diff(x, y):
    return abs(np.mod(x - y + np.pi, 2 * np.pi) - np.pi)


def random_orbit(rng):
    return dict(a=7000 + 50000 * rng.random(), e=0.01 + 0.89 * rng.random(),
                i=np.deg2rad(5 + 170 * rng.random()), raan=2 * np.pi * rng.random(),
                argp=2 * np.pi * rng.random(), nu=2 * np.pi * rng.random())


def test_round_trip_kepler_cartesian():
    rng = np.random.default_rng(1)
    for _ in range(200):
        o = random_orbit(rng)
        r, v = kep2car(o["a"], o["e"], o["i"], o["raan"], o["argp"], o["nu"], MU)
        b = car2kep(r, v, MU)
        assert b.a == pytest.approx(o["a"], rel=1e-10)
        assert b.e == pytest.approx(o["e"], abs=1e-10)
        assert b.i == pytest.approx(o["i"], abs=1e-9)
        for got, want in ((b.raan, o["raan"]), (b.argp, o["argp"]), (b.nu, o["nu"])):
            assert angle_diff(got, want) < 1e-8


def test_vectorized_kep2car_matches_single_calls():
    nu = np.linspace(0, 2 * np.pi, 7)
    r, v = kep2car(20000, 0.5, 0.3, 1.0, 2.0, nu, MU)
    for k, n in enumerate(nu):
        rk, vk = kep2car(20000, 0.5, 0.3, 1.0, 2.0, n, MU)
        np.testing.assert_allclose(r[k], rk, rtol=1e-14)
        np.testing.assert_allclose(v[k], vk, rtol=1e-14)


def test_time_of_flight_against_propagation():
    a, e, i, raan, argp, nu0, nu1 = 30000, 0.6, 0.5, 1, 2, 0.4, 4.0
    dt = time_of_flight(a, e, nu0, nu1, MU)
    r0, v0 = kep2car(a, e, i, raan, argp, nu0, MU)
    sol = solve_ivp(lambda t, y: np.r_[y[3:], -MU * y[:3] / np.linalg.norm(y[:3])**3],
                    (0, dt), np.r_[r0, v0], method="DOP853", rtol=1e-12, atol=1e-9)
    r_expected, _ = kep2car(a, e, i, raan, argp, nu1, MU)
    assert np.linalg.norm(sol.y[:3, -1] - r_expected) < 1e-3


def test_bitangent_matches_hohmann():
    r1, r2 = 7000, 42164
    dv1, dv2, dt, *_ = bitangent_transfer(r1, 0, r2, 0, "pa", 0, MU)
    hohmann = (np.sqrt(MU / r1) * (np.sqrt(2 * r2 / (r1 + r2)) - 1)
               + np.sqrt(MU / r2) * (1 - np.sqrt(2 * r1 / (r1 + r2))))
    assert dv1 + dv2 == pytest.approx(hohmann, rel=1e-12)
    assert dt == pytest.approx(np.pi * np.sqrt(((r1 + r2) / 2)**3 / MU), rel=1e-12)


def test_plane_change_keeps_position():
    rng = np.random.default_rng(3)
    for _ in range(200):
        a, e = 10000 + 30000 * rng.random(), 0.8 * rng.random()
        i_i, i_f = np.deg2rad(5 + 80 * rng.random(2))
        raan_i, raan_f, argp_i = 2 * np.pi * rng.random(3)
        dv, argp_f, nu = plane_change(a, e, i_i, raan_i, argp_i, i_f, raan_f, MU)
        r_old, v_old = kep2car(a, e, i_i, raan_i, argp_i, nu, MU)
        r_new, v_new = kep2car(a, e, i_f, raan_f, argp_f, nu, MU)
        assert np.linalg.norm(r_new - r_old) < 1e-6 * np.linalg.norm(r_old)
        assert np.linalg.norm(v_new - v_old) == pytest.approx(dv, rel=1e-9)


def test_periapsis_change_keeps_position():
    rng = np.random.default_rng(4)
    for _ in range(100):
        a, e = 10000 + 30000 * rng.random(), 0.05 + 0.85 * rng.random()
        i, raan, argp_i, argp_f = np.deg2rad(5 + 80 * rng.random()), *(2 * np.pi * rng.random(3))
        dv, nu_i, nu_f = periapsis_change(a, e, argp_i, argp_f, MU)
        for j in range(2):
            r_old, v_old = kep2car(a, e, i, raan, argp_i, nu_i[j], MU)
            r_new, v_new = kep2car(a, e, i, raan, argp_f, nu_f[j], MU)
            assert np.linalg.norm(r_new - r_old) < 1e-6 * np.linalg.norm(r_old)
            assert np.linalg.norm(v_new - v_old) == pytest.approx(dv, rel=1e-9)
