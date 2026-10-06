"""Run the three scenarios in sequence and print the mission budget (see runMission.m)."""

from types import SimpleNamespace as NS

import numpy as np

from .core import ecliptic_to_equatorial
from .scenario1 import run_scenario1
from .scenario2 import run_scenario2
from .scenario3 import run_scenario3


def run_mission(cfg):
    """Return a namespace with sc1, sc2, sc3 and the mission budget."""
    # The heliocentric transfer fixes the departure asymptote, which defines
    # the parking orbit that the geocentric transfer has to reach.
    sc2 = run_scenario2(cfg)
    v_inf = ecliptic_to_equatorial(sc2.best.v_t1 - sc2.best.v_earth, cfg.const.obliquity)
    sc1 = run_scenario1(cfg, v_inf)
    sc3 = run_scenario3(cfg, sc1, sc2)

    m = NS(dv_geocentric=sc1.best.dv_total)
    m.dv_coplanar = m.dv_geocentric + sc3.dv_coplanar
    m.dv_noncoplanar = m.dv_geocentric + sc3.dv_noncoplanar
    m.duration_days = (sc1.best.time / 86400 + sc3.escape.time_in_soi / 86400
                       + sc2.best.tof_days + sc3.capture.time_in_soi / 86400)
    res = NS(sc1=sc1, sc2=sc2, sc3=sc3, mission=m)

    if cfg.verbose:
        print_summary(cfg, res)
    if cfg.plot.enable:
        from . import plotting                       # matplotlib only when needed
        plotting.plot_all(cfg, res)
    return res


def print_summary(cfg, res):
    r1, r2, r3, m = res.sc1, res.sc2, res.sc3, res.mission
    line = "=" * 72
    deg = np.rad2deg

    print(f"\n{line}\nSCENARIO 1 - Geocentric transfer: GTO -> parking orbit\n{line}")
    if r1.launch_sweep is not None:
        print(f"Launch RAAN sweep: geocentric dv from {r1.launch_sweep.dv.min():.3f} to "
              f"{r1.launch_sweep.dv.max():.3f} km/s -> launch RAAN = {deg(r1.initial.raan):.0f} deg")
    f = r1.final
    print(f"Parking orbit ({cfg.sc1.parking.mode}): a = {f.a:.3f} km, e = {f.e:.6f}, i = {deg(f.i):.4f} deg, "
          f"RAAN = {deg(f.raan):.4f} deg, argp = {deg(f.argp):.4f} deg")
    print(f"{'Strategy':<10}{'dv [km/s]':>11}{'time [h]':>10}{'option 1':>10}{'option 2':>10}")
    for c in r1.cases:
        t = c.dt_total / 3600
        print(f"{c.name:<10}{c.dv_total:11.4f}{np.min(t):10.2f}{t[0]:10.2f}{t[1]:10.2f}")
    print(f"Selected ({cfg.sc1.selection}): {r1.best.name}, dv = {r1.best.dv_total:.4f} km/s, "
          f"time = {r1.best.time / 3600:.2f} h")

    print(f"\n{line}\nSCENARIO 2 - Heliocentric transfer: Earth -> {cfg.sc2.target.name}\n{line}")
    rec = r2.grid_records
    print(f"Valid grid points: {rec['valid'].sum()} / {rec['valid'].size}  (long arcs allowed: {cfg.sc2.allow_long_arc})")
    print(f"{'Method':<24}{'dv':>9}{'dv1':>9}{'dv2':>9}{'TOF [d]':>10}{'arc [deg]':>11}")
    for s in r2.comparison:
        print(f"{s.method:<24}{s.dv:9.4f}{s.dv1:9.4f}{s.dv2:9.4f}{s.tof_days:10.2f}{s.arc_deg:11.2f}")
    print(f"Differential evolution runs (seed, J): "
          + ", ".join(f"({int(sd)}, {J:.4f})" for sd, J in r2.de_runs))
    print(f"Multi-start: {len(r2.local_minima)} distinct valid local minima")
    b = r2.best
    print(f"Best ({b.method}): dv = {b.dv:.4f} km/s ({b.dv1:.4f} + {b.dv2:.4f}), TOF = {b.tof_days:.1f} days, "
          f"arc = {b.arc_deg:.1f} deg")

    print(f"\n{line}\nSCENARIO 3 - Patched conics\n{line}")
    print(f"SOI radius: Earth {r3.soi_earth:.0f} km, {cfg.sc2.target.name} {r3.soi_target:.3f} km")
    print(f"Escape (coplanar): v_inf = {r3.v_inf_dep:.4f} km/s, dv = {r3.escape.dv:.4f} km/s, e = {r3.escape.e:.4f}")
    print(f"Escape (non-coplanar, burn at nu = {deg(r3.escape_3d.nu_burn_on_parking):.1f} deg): "
          f"dv = {r3.escape_3d.dv:.4f} km/s")
    print(f"Capture: v_inf = {r3.v_inf_arr_mag:.4f} km/s, dv (h = {r3.capture.altitude_km} km) = "
          f"{r3.capture.dv:.4f} km/s")

    print(f"\n{line}\nMISSION BUDGET\n{line}")
    print(f"Geocentric transfer ({r1.best.name})        {m.dv_geocentric:8.4f} km/s")
    print(f"TOTAL (coplanar, idealized)          {m.dv_coplanar:8.4f} km/s")
    print(f"TOTAL (non-coplanar, realistic)      {m.dv_noncoplanar:8.4f} km/s")
    print(f"Mission duration                     {m.duration_days:8.1f} days")
