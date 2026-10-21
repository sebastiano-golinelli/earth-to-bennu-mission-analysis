"""Earth -> Bennu Mission Explorer: interactive web app on top of the bennu_mission package.

Run locally from the repository root:
    streamlit run app/streamlit_app.py
"""

import sys
from pathlib import Path
from typing import NamedTuple

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from bennu_mission import mission_config, run_mission  # noqa: E402
from bennu_mission.core import TWO_PI, Orbit, conic_points  # noqa: E402
from bennu_mission.scenario2 import evaluate_transfer, evaluate_transfers, objective_and_gradient  # noqa: E402

D2R = np.pi / 180
REPO_URL = "https://github.com/sebastiano-golinelli/earth-to-bennu-mission-analysis"
C_EARTH, C_TARGET, C_TRANSFER = "#0072BD", "#9a9a9a", "#e8432e"
PLOT_CONFIG = {"scrollZoom": False, "displaylogo": False}   # page scroll must not zoom the plots

DEFAULT = mission_config()
AU = DEFAULT.const.AU
BENNU = DEFAULT.sc2.target.orbit


class Params(NamedTuple):
    """Everything the user can change; hashable, so results can be cached."""

    a_au: float
    e: float
    i_deg: float
    raan_deg: float
    argp_deg: float
    min_perihelion_au: float
    allow_long_arc: bool
    tof_penalty: float            # [m/s per day of flight]


def make_cfg(p: Params):
    cfg = mission_config()
    cfg.verbose = False
    cfg.plot.enable = False
    s = cfg.sc2
    s.target.orbit = Orbit(a=p.a_au * AU, e=p.e, i=p.i_deg * D2R, raan=p.raan_deg * D2R, argp=p.argp_deg * D2R)
    if not _is_bennu(p):
        s.target.name = "Custom target"
    s.min_perihelion = p.min_perihelion_au * AU
    s.allow_long_arc = p.allow_long_arc
    s.weight_tof = p.tof_penalty / 1000
    # lighter global search than the reference run: same optimum for Bennu, ~5 s instead of ~25 s
    s.de.seeds = [1]
    s.multistart.points = 40
    return cfg


def _is_bennu(p: Params):
    ref = (BENNU.a / AU, BENNU.e, BENNU.i / D2R, BENNU.raan / D2R, BENNU.argp / D2R)
    return np.allclose((p.a_au, p.e, p.i_deg, p.raan_deg, p.argp_deg), ref, rtol=1e-9, atol=1e-9)


# ---- computations (cached) -----------------------------------------------------------

@st.cache_data(show_spinner="Evaluating 72,000 candidate transfers and refining the best one...")
def optimize(p: Params):
    """Vectorized grid search + SLSQP refinement (the fast part of Scenario 2)."""
    cfg = make_cfg(p)
    s = cfg.sc2
    vecs = [np.linspace(0, TWO_PI, n + 1)[:-1] for n in s.grid]
    N1, N2, N3 = np.meshgrid(*vecs, indexing="ij")
    X = np.stack([N1.ravel(), N2.ravel(), N3.ravel()], 1)
    rec = evaluate_transfers(X, cfg)
    J = np.where(rec["valid"], rec["J"], np.inf)
    if not np.isfinite(J).any():
        return None
    grid = evaluate_transfer(X[np.argmin(J)], cfg, "Grid search")
    res = minimize(lambda x: objective_and_gradient(x, cfg), grid.x, method="SLSQP", jac=True,
                   options=dict(maxiter=s.local.max_iterations, ftol=1e-14))
    best = evaluate_transfer(res.x, cfg, "Grid search + SLSQP")
    if not best.valid or best.J > grid.J:
        best = grid

    dv = np.where(rec["valid"], rec["dv"], np.inf).reshape(s.grid).min(axis=2)
    dv[np.isinf(dv)] = np.nan
    v = rec["valid"]
    return dict(grid=grid, best=best, heat=dv, vecs=[np.rad2deg(x) for x in vecs], n_valid=int(v.sum()),
                n_total=int(v.size), tof=rec["tof_days"][v], dv=rec["dv"][v], arc=rec["arc_deg"][v])


@st.cache_data(show_spinner="Running the full mission: launch orbit, heliocentric transfer, escape and capture...")
def full_mission(p: Params):
    res = run_mission(make_cfg(p))
    r1, r2, r3, m = res.sc1, res.sc2, res.sc3, res.mission
    return dict(
        launch_raan=float(np.rad2deg(r1.initial.raan)), sc1_name=r1.best.name, sc1_dv=float(r1.best.dv_total),
        sc1_hours=float(r1.best.time / 3600),
        sc1_table=[dict(Strategy=c.name, **{"Δv [km/s]": round(float(c.dv_total), 4),
                                            "Time [h]": round(float(np.min(c.dt_total) / 3600), 2)})
                   for c in sorted(r1.cases, key=lambda c: c.dv_total)],
        best=r2.best, n_minima=len(r2.local_minima),
        escape_2d=float(r3.escape.dv), escape_3d=float(r3.escape_3d.dv), v_inf_dep=r3.v_inf_dep,
        capture=float(r3.capture.dv), capture_alt=float(r3.capture.altitude_km), v_inf_arr=r3.v_inf_arr_mag,
        soi_target=float(r3.soi_target), total=float(m.dv_noncoplanar), total_2d=float(m.dv_coplanar),
        days=float(m.duration_days))


# ---- plots -----------------------------------------------------------------------------

def orbit_figure(cfg, sol, z_scale=1.0):
    nu = np.linspace(0, TWO_PI, 500)
    oE, oT = cfg.sc2.earth, cfg.sc2.target.orbit
    scale = np.array([1, 1, z_scale])
    PE = conic_points(oE.a, oE.e, oE.i, oE.raan, oE.argp, nu) / AU * scale
    PT = conic_points(oT.a, oT.e, oT.i, oT.raan, oT.argp, nu) / AU * scale
    fig = go.Figure()

    def line(P, name, color, width, dash=None):
        fig.add_scatter3d(x=P[:, 0], y=P[:, 1], z=P[:, 2], mode="lines", name=name,
                          line=dict(color=color, width=width, dash=dash), hoverinfo="name")

    def point(r, name, color, size, symbol="circle"):
        fig.add_scatter3d(x=[r[0]], y=[r[1]], z=[r[2]], mode="markers", name=name, hoverinfo="name",
                          marker=dict(color=color, size=size, symbol=symbol, line=dict(color="black", width=1)))

    line(PE, "Earth orbit", C_EARTH, 4)
    line(PT, f"{cfg.sc2.target.name} orbit", C_TARGET, 4)
    if sol.valid:
        arc = sol.nu1_t + np.linspace(0, np.deg2rad(sol.arc_deg), 400)
        PX = conic_points(sol.a_t, sol.e_t, sol.i_t, sol.raan_t, sol.argp_t, arc) / AU * scale
        line(PX, "Transfer", C_TRANSFER, 7)
    point([0, 0, 0], "Sun", "#ffcc00", 9)
    point(sol.r1 / AU * scale, "Departure (Earth)", C_EARTH, 6)
    point(sol.r2 / AU * scale, "Arrival", C_TARGET, 6)
    zlabel = "z [AU]" if z_scale == 1 else f"z [AU] (× {z_scale:g})"
    fig.update_layout(height=540, margin=dict(l=0, r=0, t=10, b=0), legend=dict(orientation="h", y=1.02),
                      scene=dict(aspectmode="data", xaxis_title="x [AU]", yaxis_title="y [AU]", zaxis_title=zlabel,
                                 camera=dict(eye=dict(x=0.9, y=-1.3, z=0.8))))
    return fig


def heatmap_figure(opt):
    best = opt["best"]
    fig = go.Figure(go.Heatmap(x=opt["vecs"][1], y=opt["vecs"][0], z=opt["heat"], colorscale="Viridis",
                               zmin=np.nanmin(opt["heat"]), zmax=np.nanmin(opt["heat"]) + 10,
                               colorbar=dict(title="Δv [km/s]"),
                               hovertemplate="arrival %{x:.0f}°<br>departure %{y:.0f}°<br>Δv %{z:.3f} km/s<extra></extra>"))
    fig.add_scatter(x=[np.rad2deg(best.x[1]) % 360], y=[np.rad2deg(best.x[0]) % 360], mode="markers",
                    marker=dict(symbol="star", size=18, color="red", line=dict(color="white", width=1)),
                    name="Best transfer", hoverinfo="name")
    fig.update_layout(height=430, margin=dict(l=0, r=0, t=30, b=0), showlegend=False,
                      title="Minimum Δv for each departure / arrival position",
                      xaxis_title="Arrival true anomaly on the target orbit [deg]",
                      yaxis_title="Departure true anomaly on Earth's orbit [deg]")
    return fig


def tradeoff_figure(opt):
    best = opt["best"]
    keep = opt["dv"] <= best.dv + 5
    tof, dv, arc = opt["tof"][keep], opt["dv"][keep], opt["arc"][keep]
    if tof.size > 15000:                                   # keep the browser responsive
        idx = np.random.default_rng(0).choice(tof.size, 15000, replace=False)
        tof, dv, arc = tof[idx], dv[idx], arc[idx]
    fig = go.Figure(go.Scattergl(x=tof, y=dv, mode="markers", name="Grid transfers",
                                 marker=dict(size=4, color=arc, colorscale="Viridis",
                                             colorbar=dict(title="Transfer<br>angle [deg]")),
                                 hovertemplate="TOF %{x:.0f} d<br>Δv %{y:.3f} km/s<extra></extra>"))
    fig.add_scatter(x=[best.tof_days], y=[best.dv], mode="markers", name="Best transfer",
                    marker=dict(symbol="star", size=18, color="yellow", line=dict(color="black", width=1)))
    fig.update_layout(height=430, margin=dict(l=0, r=0, t=30, b=0), showlegend=False,
                      title="Cost versus time of flight", xaxis_title="Time of flight [days]",
                      yaxis_title="Δv departure + arrival [km/s]")
    return fig


def transfer_metrics(sol):
    c = st.columns(5)
    c[0].metric("Total Δv", f"{sol.dv:.3f} km/s")
    c[1].metric("Departure Δv", f"{sol.dv1:.3f} km/s")
    c[2].metric("Arrival Δv", f"{sol.dv2:.3f} km/s")
    c[3].metric("Time of flight", f"{sol.tof_days:.0f} days")
    c[4].metric("Transfer angle", f"{sol.arc_deg:.0f}°")


def why_invalid(sol, cfg):
    if not np.isfinite(sol.e_t):
        return "Departure and arrival points are aligned with the Sun: the transfer plane is undefined."
    if sol.e_t >= 1:
        return "Joining these two points with this periapsis direction needs a hyperbola, not an ellipse."
    if sol.e_t < 0:
        return "No ellipse with this periapsis direction passes through both points."
    if np.isfinite(sol.rp_t) and sol.rp_t <= cfg.sc2.min_perihelion:
        return (f"The transfer would pass {sol.rp_t / AU:.2f} AU from the Sun, below the "
                f"{cfg.sc2.min_perihelion / AU:.2f} AU limit.")
    return "This combination does not give a valid transfer."


# ---- page ------------------------------------------------------------------------------

st.set_page_config(page_title="Earth → Bennu Mission Explorer", page_icon="🛰️", layout="wide")
st.title("🛰️ Earth → Bennu Mission Explorer")
st.caption("Preliminary design of a mission to the near-Earth asteroid 101955 Bennu: interactive version of "
           f"an open-source [MATLAB + Python mission analysis]({REPO_URL}).")

with st.sidebar:
    st.header("Mission settings")
    with st.expander("Target orbit (default: Bennu, JPL SBDB)"):
        if st.button("Reset to Bennu"):
            for k in ("a", "e", "i", "raan", "argp"):
                st.session_state.pop(f"tgt_{k}", None)
        a_au = st.number_input("Semi-major axis [AU]", 0.6, 3.0, BENNU.a / AU, 0.01, format="%.4f", key="tgt_a")
        ecc = st.number_input("Eccentricity [-]", 0.0, 0.9, BENNU.e, 0.01, format="%.4f", key="tgt_e")
        inc = st.number_input("Inclination [deg]", 0.0, 60.0, BENNU.i / D2R, 0.5, format="%.3f", key="tgt_i")
        raan = st.number_input("RAAN [deg]", 0.0, 360.0, BENNU.raan / D2R, 1.0, format="%.3f", key="tgt_raan")
        argp = st.number_input("Argument of perihelion [deg]", 0.0, 360.0, BENNU.argp / D2R, 1.0, format="%.3f",
                               key="tgt_argp")
    min_peri = st.slider("Minimum distance from the Sun [AU]", 0.1, 0.9, 0.3, 0.05,
                         help="Transfers that pass closer to the Sun are rejected (thermal limit).")
    long_arc = st.checkbox("Allow transfers longer than half an orbit", True)
    tof_penalty = st.slider("Penalty on flight time [m/s per day]", 0.0, 20.0, 0.0, 0.5,
                            help="0 = cheapest transfer. Higher values trade Δv for a faster trip.")
    z_scale = st.select_slider("Out-of-plane exaggeration in 3D plots", [1, 5, 10, 20], 1)
    st.divider()
    st.markdown(f"Code and validation: [GitHub]({REPO_URL})")

params = Params(a_au, ecc, inc, raan, argp, min_peri, long_arc, tof_penalty)
cfg = make_cfg(params)

tab_opt, tab_manual, tab_full, tab_about = st.tabs(
    ["Optimal transfer", "Design it yourself", "Full mission budget", "About the model"])

with tab_opt:
    opt = optimize(params)
    if opt is None:
        st.error("No valid transfer exists with these settings. Try a lower Sun-distance limit or other target "
                 "elements.")
    else:
        best = opt["best"]
        st.subheader(f"Cheapest direct transfer to {cfg.sc2.target.name}")
        transfer_metrics(best)
        st.plotly_chart(orbit_figure(cfg, best, z_scale), width="stretch", config=PLOT_CONFIG)
        left, right = st.columns(2)
        left.plotly_chart(heatmap_figure(opt), width="stretch", config=PLOT_CONFIG)
        right.plotly_chart(tradeoff_figure(opt), width="stretch", config=PLOT_CONFIG)
        st.caption(f"{opt['n_valid']:,} of {opt['n_total']:,} grid points give a valid elliptic transfer. "
                   "The best grid point is refined with a gradient-based optimizer (SLSQP); the full-mission tab "
                   "also runs a global search (differential evolution and multi-start).")

with tab_manual:
    st.subheader("Pick the departure and arrival points yourself")
    st.write("A two-impulse transfer is defined by where the spacecraft leaves Earth's orbit, where it meets the "
             "target, and the orientation of the transfer ellipse. Move the sliders and see what it costs.")

    def load_optimum():
        if opt is not None:
            x = np.rad2deg(np.mod(opt["best"].x, TWO_PI))
            st.session_state.update(nu1=float(round(x[0])), nu2=float(round(x[1])), w=float(round(x[2])))

    if "nu1" not in st.session_state:
        load_optimum()
    st.button("Start from the optimum", on_click=load_optimum)
    c1, c2, c3 = st.columns(3)
    nu1 = c1.slider("Departure position on Earth's orbit [deg]", 0.0, 359.0, step=1.0, format="%.0f°", key="nu1")
    nu2 = c2.slider("Arrival position on the target orbit [deg]", 0.0, 359.0, step=1.0, format="%.0f°", key="nu2")
    w = c3.slider("Orientation of the transfer ellipse [deg]", 0.0, 359.0, step=1.0, format="%.0f°", key="w")
    sol = evaluate_transfer(np.deg2rad([nu1, nu2, w]), cfg, "Manual")
    if sol.valid:
        transfer_metrics(sol)
        if opt is not None:
            extra = sol.dv - opt["best"].dv
            st.info(f"This transfer costs {1000 * extra:.0f} m/s more than the optimum." if extra >= 5e-4
                    else "Within 1 m/s of the optimum.")
    else:
        st.warning(why_invalid(sol, cfg))
    st.plotly_chart(orbit_figure(cfg, sol, z_scale), width="stretch", config=PLOT_CONFIG)

with tab_full:
    st.subheader("From the launch orbit to orbit around the asteroid")
    st.write("Chains the three phases of the mission: transfer from the launch orbit (GTO) to a parking orbit "
             "designed for the escape, heliocentric transfer, and capture into a 1 km orbit at the target.")
    if st.button("Run the full mission (about 5 seconds)", type="primary") or st.session_state.get("ran_full"):
        st.session_state.ran_full = True
        try:
            f = full_mission(params)
        except RuntimeError as err:
            st.error(f"The mission cannot be closed with these settings: {err}")
        else:
            c = st.columns(3)
            c[0].metric("Total Δv", f"{f['total']:.3f} km/s")
            c[1].metric("Mission duration", f"{f['days']:.0f} days")
            c[2].metric("Launch window (RAAN)", f"{f['launch_raan']:.0f}°")
            b = f["best"]
            st.table([
                {"Phase": "GTO → parking orbit", "Δv [km/s]": f"{f['sc1_dv']:.4f}",
                 "Notes": f"Strategy {f['sc1_name']}, {f['sc1_hours']:.1f} h"},
                {"Phase": "Escape from Earth", "Δv [km/s]": f"{f['escape_3d']:.4f}",
                 "Notes": f"v∞ = {f['v_inf_dep']:.3f} km/s (coplanar estimate {f['escape_2d']:.4f})"},
                {"Phase": "Heliocentric cruise", "Δv [km/s]": "—",
                 "Notes": f"{b.tof_days:.0f} days, transfer angle {b.arc_deg:.0f}°"},
                {"Phase": "Capture at the target", "Δv [km/s]": f"{f['capture']:.4f}",
                 "Notes": f"v∞ = {f['v_inf_arr']:.3f} km/s, {f['capture_alt']:g} km altitude, "
                          f"SOI {f['soi_target']:.2f} km"},
                {"Phase": "Total", "Δv [km/s]": f"{f['total']:.4f}", "Notes": f"{f['days']:.0f} days"},
            ])
            with st.expander("Geocentric strategies compared (S = shape, P = plane, W = periapsis)"):
                st.dataframe(f["sc1_table"], hide_index=True, width="stretch")
            st.caption(f"Global search found {f['n_minima']} distinct local minima of the heliocentric problem.")

with tab_about:
    st.markdown(f"""
**What this is.** A preliminary (patched-conics) mission analysis from a geostationary transfer orbit to the
near-Earth asteroid Bennu, with orbital data from the JPL Small-Body Database. The model and every number on this
page come from the open-source [`bennu_mission`]({REPO_URL}) package.

**How it is validated.** The same analysis exists in MATLAB and in Python. The Python results match the MATLAB
ones to 1e-14 relative error, and both versions are checked by automated tests on every commit (GitHub Actions):
analytical cases, numerical propagation of the transfer, and regression against reference values.

**Model limits.** Two-body arcs, impulsive burns, no launch date (orbits are fixed, so the result is the best
geometry, not a calendar launch window), no planetary flybys. Good for a first sizing, not for flight.

**Tech.** Python, NumPy, SciPy (SLSQP, differential evolution), Plotly, Streamlit. The 72,000-point grid search is
fully vectorized and runs in well under a second.

Built by Sebastiano Enrico Golinelli, aerospace engineering student at Politecnico di Milano.
""")
