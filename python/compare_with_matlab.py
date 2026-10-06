"""Print how closely the Python port reproduces the MATLAB results.

Deterministic quantities are computed from the MATLAB optimal transfer, so
that they depend only on the ported formulas; the optimizers are compared on
their final optimum. Run from the python/ folder: python compare_with_matlab.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "tests"))
from conftest import ROOT  # noqa: E402,F401  (adds the package to the path)
import conftest  # noqa: E402

from bennu_mission import mission_config, run_mission  # noqa: E402
from bennu_mission.scenario2 import evaluate_transfers  # noqa: E402
from test_matlab_parity import grid_points  # noqa: E402


def max_rel(actual, expected):
    a, e = np.asarray(actual, float).ravel(), np.asarray(expected, float).ravel()
    ok = np.isfinite(e)
    assert np.array_equal(ok, np.isfinite(a)), "NaN pattern differs"
    return float(np.max(np.abs(a[ok] - e[ok]) / np.abs(e[ok]))), int(ok.sum())


def main():
    cfg = mission_config()
    cfg.plot.enable = False
    cfg.verbose = False
    ref = conftest.ref.__wrapped__()
    det = conftest.from_matlab_optimum.__wrapped__(cfg, ref)
    full = run_mission(cfg)
    r1, r3, s1, s2, s3 = det.sc1, det.sc3, ref["sc1"], ref["sc2"], ref["sc3"]

    rows = [
        ("Launch window sweep (dv)", *max_rel(r1.launch_sweep.dv, s1["launchSweepDv"])),
        ("16 geocentric strategies (dv)", *max_rel([c.dv_total for c in r1.cases], s1["dv"])),
        ("16 geocentric strategies (time)", *max_rel([np.min(c.dt_total) / 3600 for c in r1.cases], s1["timeH"])),
        ("Parking orbit elements", *max_rel([r1.final.a, r1.final.e, r1.final.i, r1.final.raan, r1.final.argp],
                                            s1["parking"])),
        ("Grid search sample (dv)", *max_rel(evaluate_transfers(grid_points(cfg, s2["gridSampleIndex"]), cfg)["dv"],
                                             s2["gridSampleDv"])),
        ("Transfer at the optimum (dv, TOF)", *max_rel([det.sc2.best.dv, det.sc2.best.tof_days],
                                                       [s2["bestDv"], s2["bestTofDays"]])),
        ("Escape, coplanar (dv)", *max_rel(r3.escape.dv, s3["escapeDv"])),
        ("Escape, 3D burn sweep (dv)", *max_rel(r3.burn_sweep.dv, s3["burnSweepDv"])),
        ("Capture table (dv)", *max_rel([row["dv_circular"] for row in r3.capture_table], s3["captureDvCircular"])),
    ]
    print(f"\n{'Deterministic quantity':<36}{'values':>8}{'max rel. difference':>22}")
    for name, err, n in rows:
        print(f"{name:<36}{n:>8}{err:>22.1e}")

    print(f"\n{'Optimizer (Python vs MATLAB)':<36}{'Python dv':>14}{'MATLAB dv':>14}{'difference':>12}")
    pairs = [("SLSQP vs fmincon", full.sc2.local.dv, s2["fminconDv"]),
             ("DE + SLSQP vs ga + fmincon", full.sc2.de_refined.dv, s2["bestDv"]),
             ("Multi-start vs MultiStart", full.sc2.multistart.dv, s2["multiStartDv"])]
    for name, py, ml in pairs:
        print(f"{name:<36}{py:14.7f}{ml:14.7f}{py - ml:12.1e}")
    print(f"\nMission total (realistic): Python {full.mission.dv_noncoplanar:.6f} km/s, "
          f"MATLAB {ref['mission']['dvNonCoplanar']:.6f} km/s")


if __name__ == "__main__":
    main()
