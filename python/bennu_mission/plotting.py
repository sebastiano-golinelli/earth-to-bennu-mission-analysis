"""Matplotlib figures (counterpart of src/plotting in MATLAB)."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .core import conic_points

C_EARTH, C_BENNU, C_PARK = "#0072BD", "#777777", "#D95319"


def plot_all(cfg, res):
    figs = {}
    figs.update(plot_scenario1(cfg, res.sc1))
    figs.update(plot_scenario2(cfg, res.sc2))
    figs.update(plot_scenario3(cfg, res.sc3))
    if cfg.plot.export:
        folder = Path(cfg.plot.folder)
        folder.mkdir(parents=True, exist_ok=True)
        for name, fig in figs.items():
            fig.savefig(folder / f"{name}.png", dpi=150, bbox_inches="tight")
    else:
        plt.show()
    return figs


def _earth(ax, radius):
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 40), np.linspace(0, np.pi, 20))
    ax.plot_surface(radius * np.cos(u) * np.sin(v), radius * np.sin(u) * np.sin(v), radius * np.cos(v),
                    color="#6fa8dc", alpha=0.5, linewidth=0)


def _equal_3d(ax, P):
    center = (P.max(0) + P.min(0)) / 2
    half = (P.max(0) - P.min(0)).max() / 2
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(center[2] - half, center[2] + half)
    ax.set_box_aspect((1, 1, 1))


def _orbit(o, n=600):
    return conic_points(o.a, o.e, o.i, o.raan, o.argp, np.linspace(0, 2 * np.pi, n))


def plot_scenario1(cfg, r1):
    figs = {}
    o_i, o_f = r1.initial, r1.final

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(projection="3d")
    _earth(ax, cfg.const.r_earth)
    P_i, P_f = _orbit(o_i), _orbit(o_f)
    ax.plot(*P_i.T, color=C_EARTH, label="Launch orbit (GTO)")
    ax.plot(*P_f.T, color=C_PARK, label="Parking orbit")
    burn = conic_points(o_f.a, o_f.e, o_f.i, o_f.raan, o_f.argp, 0.0)
    ax.scatter(*burn, color="red", marker="*", s=150, label="Escape burn point", depthshade=False)
    _equal_3d(ax, np.vstack([P_i, P_f]))
    ax.set(xlabel="x [km]", ylabel="y [km]", zlabel="z [km]", title="Launch and parking orbits (equatorial frame)")
    ax.legend(loc="upper left")
    figs["sc1_orbits"] = fig

    if r1.launch_sweep is not None:
        s = r1.launch_sweep
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(s.raan_deg, s.dv, lw=1.6, label="Best strategy")
        ax.plot(np.rad2deg(o_i.raan), r1.best.dv_total, "r*", ms=14,
                label=f"Selected: {np.rad2deg(o_i.raan):.0f} deg, {r1.best.dv_total:.3f} km/s")
        ax.set(xlim=(0, 360), xticks=range(0, 361, 60), xlabel="RAAN of the launch orbit [deg] (set by the launch time)",
               ylabel=r"Geocentric $\Delta v$ [km/s]", title="Launch window: cheapest strategy for each launch RAAN")
        ax.grid(alpha=0.4)
        ax.legend(loc="upper center")
        figs["sc1_launch_window"] = fig

    strategies, types = cfg.sc1.strategies, cfg.sc1.types
    dv = np.array([c.dv_total for c in r1.cases]).reshape(len(strategies), len(types))
    time = np.array([np.min(c.dt_total) / 3600 for c in r1.cases]).reshape(len(strategies), len(types))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, data, fmt, title in ((axes[0], dv, "{:.3f}", r"$\Delta v$ [km/s]"),
                                 (axes[1], time, "{:.1f}", "Transfer time [h]")):
        ax.imshow(data, cmap="Blues")
        for (k, j), val in np.ndenumerate(data):
            ax.text(j, k, fmt.format(val), ha="center", va="center",
                    color="white" if val > np.nanmean(data) else "black")
        ax.set(xticks=range(len(types)), xticklabels=types, yticks=range(len(strategies)),
               yticklabels=strategies, title=title, xlabel="Bitangent (departure, arrival apsis)")
    axes[0].set_ylabel("Strategy")
    fig.suptitle(f"Geocentric strategies (S = shape, P = plane, W = periapsis) - selected: {r1.best.name}")
    fig.tight_layout()
    figs["sc1_strategies"] = fig
    return figs


def plot_scenario2(cfg, r2):
    figs = {}
    AU = cfg.const.AU
    rec, best = r2.grid_records, r2.best

    n1, n2, n3 = cfg.sc2.grid
    dv = np.where(rec["valid"], rec["dv"], np.nan).reshape(n1, n2, n3)
    with np.errstate(all="ignore"):
        dv_min = np.nanmin(np.where(np.isnan(dv), np.inf, dv), axis=2)
    dv_min[np.isinf(dv_min)] = np.nan
    fig, ax = plt.subplots(figsize=(8, 5.5))
    h1, h2 = 180 / n1, 180 / n2                 # half cell: pixels centred on the grid values
    im = ax.imshow(dv_min, origin="lower", aspect="auto", cmap="viridis",
                   extent=[-h2, 360 - h2, -h1, 360 - h1], vmin=np.nanmin(dv_min), vmax=np.nanmin(dv_min) + 10)
    ax.plot(np.rad2deg(best.x[1]), np.rad2deg(best.x[0]), "r*", ms=14)
    fig.colorbar(im, ax=ax, label=r"min $\Delta v$ over $\omega_T$ [km/s]")
    ax.set(xlabel="Arrival true anomaly on Bennu's orbit [deg]", ylabel="Departure true anomaly on Earth's orbit [deg]",
           title=r"Grid search: minimum $\Delta v$ over $\omega_T$")
    figs["sc2_grid_search"] = fig

    v = rec["valid"]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    sc = ax.scatter(rec["tof_days"][v], rec["dv"][v], s=3, c=rec["arc_deg"][v], cmap="viridis", label="Grid points")
    fig.colorbar(sc, ax=ax, label="Transfer angle [deg]")
    lm = r2.local_minima
    ax.plot([m["tof_days"] for m in lm], [m["dv"] for m in lm], "ko", mfc="none", ms=7, label="Multi-start local minima")
    ax.plot(best.tof_days, best.dv, "*", ms=16, mfc="yellow", mec="k", label="Best transfer")
    # x range from the points inside the y window only (very eccentric transfers
    # with huge times of flight are valid but cost far more)
    shown = v & (rec["dv"] <= best.dv + 5)
    ax.set(ylim=(0.95 * best.dv, best.dv + 5), xlim=(0, 1.05 * rec["tof_days"][shown].max()),
           xlabel="Time of flight [days]", ylabel=r"$\Delta v_1 + \Delta v_2$ [km/s]",
           title="Heliocentric transfers: cost versus time of flight")
    ax.grid(alpha=0.4)
    ax.legend(loc="upper right")
    figs["sc2_pareto"] = fig

    oE, oT = cfg.sc2.earth, cfg.sc2.target.orbit
    PE, PB = _orbit(oE) / AU, _orbit(oT) / AU
    PT = conic_points(best.a_t, best.e_t, best.i_t, best.raan_t, best.argp_t,
                      best.nu1_t + np.linspace(0, np.deg2rad(best.arc_deg), 400)) / AU
    r1, r2pos = best.r1 / AU, best.r2 / AU
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2))
    ax1.plot(PE[:, 0], PE[:, 1], color=C_EARTH, label="Earth orbit")
    ax1.plot(PB[:, 0], PB[:, 1], color=C_BENNU, label="Bennu orbit")
    ax1.plot(PT[:, 0], PT[:, 1], "r", lw=2.4, label="Transfer")
    ax1.plot(0, 0, "o", ms=12, mfc="#ffcc00", mec="k", label="Sun")
    ax1.plot(*r1[:2], "o", ms=8, mfc=C_EARTH, mec="k", label="Departure")
    ax1.plot(*r2pos[:2], "o", ms=8, mfc=C_BENNU, mec="k", label="Arrival")
    ax1.set(aspect="equal", xlabel="x [AU]", ylabel="y [AU]", title="Ecliptic plane, top view")
    ax1.grid(alpha=0.4)
    for P, style in ((PE, dict(color=C_EARTH)), (PB, dict(color=C_BENNU)), (PT, dict(color="r", lw=2.4))):
        lon, z = _longitude_height(P)
        ax2.plot(lon, z, **style)
    for p, col in ((r1, C_EARTH), (r2pos, C_BENNU)):
        lon, z = _longitude_height(p[None, :])
        ax2.plot(lon, z, "o", ms=8, mfc=col, mec="k")
    ax2.set(xlim=(0, 360), xticks=range(0, 361, 60), xlabel="Ecliptic longitude [deg]",
            ylabel="Height above the ecliptic [AU]", title="Out-of-plane motion: arrival at Bennu's node")
    ax2.grid(alpha=0.4)
    fig.legend(*ax1.get_legend_handles_labels(), loc="lower center", ncol=6)
    fig.suptitle(rf"Best transfer: $\Delta v$ = {best.dv:.3f} km/s, TOF = {best.tof_days:.0f} days, "
                 f"transfer angle {best.arc_deg:.0f} deg")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    figs["sc2_best_transfer"] = fig

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    hist = r2.local_history
    it = np.arange(hist.size)
    bad = hist >= 1e8                       # trial points in the invalid region (penalty)
    axes[0].plot(it[~bad], hist[~bad], "o-", label="Valid iterate")
    top = hist[~bad].max() + 0.2 * np.ptp(hist[~bad])
    axes[0].plot(it[bad], np.full(bad.sum(), top), "rx", ms=9, mew=2, label="Invalid iterate (penalty)")
    axes[0].set(xlabel="Iteration", ylabel=r"$\Delta v$ [km/s]", title="SLSQP from the best grid point")
    axes[0].legend(loc="center right")
    axes[1].plot(np.arange(1, r2.de_history.size + 1), r2.de_history)
    axes[1].set(xlabel="Generation", ylabel=r"Best $\Delta v$ [km/s]", title="Differential evolution, best seed",
                ylim=(best.dv - 0.05, min(r2.de_history.max(), best.dv + 2)))
    axes[2].bar([f"seed {int(s)}" for s in r2.de_runs[:, 0]], r2.de_runs[:, 1])
    axes[2].axhline(best.dv, color="r", ls="--", label="Global optimum")
    axes[2].set(ylim=(best.dv - 0.1, r2.de_runs[:, 1].max() + 0.1), ylabel=r"Final $\Delta v$ [km/s]",
                title="Differential evolution: result of each seed")
    axes[2].legend(loc="upper left")
    for ax in axes:
        ax.grid(alpha=0.4)
    fig.tight_layout()
    figs["sc2_optimizers"] = fig
    return figs


def _longitude_height(P):
    """Ecliptic longitude [deg] and height of (N, 3) points, NaN where the longitude wraps."""
    lon = np.mod(np.rad2deg(np.arctan2(P[:, 1], P[:, 0])), 360)
    z = P[:, 2].astype(float)
    jumps = np.flatnonzero(np.abs(np.diff(lon)) > 180) + 1
    return np.insert(lon, jumps, np.nan), np.insert(z, jumps, np.nan)


def plot_scenario3(cfg, r3):
    figs = {}
    park, H = r3.parking, r3.escape_3d

    fig = plt.figure(figsize=(8, 6.5))
    ax = fig.add_subplot(projection="3d")
    _earth(ax, cfg.const.r_earth)
    P = _orbit(park)
    ax.plot(*P.T, color=C_PARK, label="Parking orbit")
    r_max = 2.5 * park.ra
    Q = conic_points(H.a, H.e, H.i, H.raan, H.argp, np.linspace(H.nu_burn, 0.98 * H.nu_inf, 2000))
    Q = Q[np.linalg.norm(Q, axis=1) <= r_max]
    ax.plot(*Q.T, "r", lw=2, label="Escape hyperbola")
    ax.scatter(*H.r_burn, color="yellow", edgecolor="k", marker="*", s=200, label="Escape burn", depthshade=False)
    d = r3.v_inf_dep_equatorial / np.linalg.norm(r3.v_inf_dep_equatorial) * r_max
    ax.quiver(0, 0, 0, *d, color="m", lw=1.6, arrow_length_ratio=0.08, label=r"Departure $v_\infty$ direction")
    _equal_3d(ax, np.vstack([P, Q, d[None, :], np.zeros((1, 3))]))
    ax.set(xlabel="x [km]", ylabel="y [km]", zlabel="z [km]",
           title=rf"Escape: $\Delta v$ = {H.dv:.4f} km/s, $v_\infty$ = {r3.v_inf_dep:.3f} km/s")
    ax.view_init(25, 35)
    ax.legend(loc="upper left")
    figs["sc3_escape"] = fig

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.semilogy(np.rad2deg(r3.burn_sweep.nu), r3.burn_sweep.dv, lw=1.6, label="3D escape from this point")
    ax.axhline(r3.escape.dv, color="k", ls="--", label="Coplanar estimate (tangential burn at periapsis)")
    ax.plot(np.rad2deg(H.nu_burn_on_parking), H.dv, "r*", ms=14, label="Selected burn point")
    ax.set(xlim=(0, 360), xticks=range(0, 361, 60), xlabel="True anomaly of the burn point on the parking orbit [deg]",
           ylabel=r"Escape $\Delta v$ [km/s]",
           title="Non-coplanar escape: cost of the burn point\n(gaps: the escape hyperbola would pass below the Earth's surface)")
    ax.grid(alpha=0.4, which="both")
    ax.legend(loc="upper left")
    figs["sc3_burn_point"] = fig

    T = r3.capture_table
    radius = np.array([r["radius_km"] for r in T])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    ax1.semilogy(radius, [r["radius_over_soi"] for r in T], "o-", label=r"$r / r_{SOI}$")
    ax1.semilogy(radius, [r["srp_ratio"] for r in T], "s-", label="Solar radiation pressure / gravity")
    ax1.semilogy(radius, [r["tide_ratio"] for r in T], "d-", label="Solar tide / gravity")
    ax1.axhline(0.1, color="k", ls=":")
    ax1.set(xlabel="Capture orbit radius [km]", ylabel="Ratio [-]", title="Is the two-body capture model valid?")
    ax1.legend(loc="lower right")
    ax2.plot(radius, [1000 * (r3.v_inf_arr_mag - r["dv_circular"]) for r in T], "o-")
    ax2.set(xlabel="Capture orbit radius [km]", ylabel=r"$v_\infty - \Delta v_{capture}$ [m/s]",
            title="Saving due to Bennu's gravity")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.4, which="both")
    fig.suptitle(rf"Capture at Bennu (SOI radius {r3.soi_target:.2f} km, $v_\infty$ = {r3.v_inf_arr_mag:.3f} km/s)")
    fig.tight_layout()
    figs["sc3_capture"] = fig
    return figs
