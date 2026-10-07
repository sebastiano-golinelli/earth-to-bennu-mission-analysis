"""Python versus MATLAB.

Deterministic parts are compared starting from the same point (the MATLAB
optimal transfer) and must agree to round-off level. The optimizers are
different algorithms, so only their final optimum is compared.
"""

import numpy as np
import pytest

from bennu_mission.scenario2 import evaluate_transfers

RTOL = 1e-9


def close(actual, expected, rtol=RTOL):
    np.testing.assert_allclose(np.asarray(actual, float), np.asarray(expected, float), rtol=rtol, atol=0,
                               equal_nan=True)


# ---- Scenario 1 -----------------------------------------------------------------

def test_launch_window_sweep(from_matlab_optimum, ref):
    sweep = from_matlab_optimum.sc1.launch_sweep
    close(sweep.dv, ref["sc1"]["launchSweepDv"])
    assert np.rad2deg(from_matlab_optimum.sc1.initial.raan) == pytest.approx(ref["sc1"]["launchRaanDeg"])


def test_strategy_table(from_matlab_optimum, ref):
    cases = from_matlab_optimum.sc1.cases
    assert [c.name for c in cases] == list(ref["sc1"]["strategies"])
    close([c.dv_total for c in cases], ref["sc1"]["dv"])
    close([np.min(c.dt_total) / 3600 for c in cases], ref["sc1"]["timeH"])
    assert from_matlab_optimum.sc1.best.name == ref["sc1"]["selected"]


def test_parking_orbit(from_matlab_optimum, ref):
    f = from_matlab_optimum.sc1.final
    close([f.a, f.e, f.i, f.raan, f.argp], ref["sc1"]["parking"])


# ---- Scenario 2 -----------------------------------------------------------------

def grid_points(cfg, index):
    n1, n2, n3 = cfg.sc2.grid
    vec = [np.linspace(0, 2 * np.pi, n + 1)[:-1] for n in (n1, n2, n3)]
    index = np.asarray(index, int)
    return np.stack([vec[0][index // (n2 * n3)], vec[1][(index // n3) % n2], vec[2][index % n3]], 1)


def test_grid_values(cfg, ref):
    out = evaluate_transfers(grid_points(cfg, ref["sc2"]["gridSampleIndex"]), cfg)
    close(out["dv"], ref["sc2"]["gridSampleDv"])


def test_grid_search_result(full_run, ref):
    rec = full_run.sc2.grid_records
    assert int(rec["valid"].sum()) == int(ref["sc2"]["gridValidCount"])
    close(full_run.sc2.grid.x, ref["sc2"]["gridBestX"])
    close(full_run.sc2.grid.dv, ref["sc2"]["gridBestDv"])


def test_transfer_at_matlab_optimum(from_matlab_optimum, ref):
    b = from_matlab_optimum.sc2.best
    close([b.dv, b.dv1, b.dv2, b.tof_days, b.arc_deg],
          [ref["sc2"][k] for k in ("bestDv", "bestDv1", "bestDv2", "bestTofDays", "bestArcDeg")])


def test_optimizers_reach_the_matlab_optimum(full_run, ref):
    # Deterministic searches must find the MATLAB optimum. Differential
    # evolution is stochastic: on another platform round-off can send a run to
    # a different basin, so only the validity of its result is required.
    r2 = full_run.sc2
    for sol in (r2.local, r2.multistart, r2.best):
        assert sol.dv == pytest.approx(ref["sc2"]["bestDv"], abs=1e-6)
    assert r2.de_refined.valid


# ---- Scenario 3 -----------------------------------------------------------------

def test_escape_and_capture(from_matlab_optimum, ref):
    r3, s = from_matlab_optimum.sc3, ref["sc3"]
    close([r3.soi_earth, r3.soi_target], [s["soiEarth"], s["soiTarget"]])
    close([r3.escape.dv, r3.escape.e, r3.escape.time_in_soi], [s["escapeDv"], s["escapeEcc"], s["escapeTimeInSoi"]])
    close([row["dv_circular"] for row in r3.capture_table], s["captureDvCircular"])
    close([row["srp_ratio"] for row in r3.capture_table], s["captureSrpRatio"])
    close([r3.dv_coplanar, r3.dv_noncoplanar], [s["dvCoplanar"], s["dvNonCoplanar"]])


def test_non_coplanar_escape_sweep(from_matlab_optimum, ref):
    r3 = from_matlab_optimum.sc3
    close(r3.burn_sweep.dv, ref["sc3"]["burnSweepDv"], rtol=1e-8)
    close(r3.escape_3d.dv, ref["sc3"]["escape3dDv"], rtol=1e-8)
    assert np.rad2deg(r3.escape_3d.nu_burn_on_parking) == pytest.approx(ref["sc3"]["burnNuDeg"])


def test_mission_budget(full_run, ref):
    m = full_run.mission
    for key, value in (("dvCoplanar", m.dv_coplanar), ("dvNonCoplanar", m.dv_noncoplanar),
                       ("durationDays", m.duration_days)):
        assert value == pytest.approx(ref["mission"][key], rel=1e-5)
