"""Earth -> Bennu Mission Explorer: interactive web app on top of the bennu_mission package.

Run locally from the repository root (so that .streamlit/config.toml is found):
    streamlit run app/streamlit_app.py
"""

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from explorer import AU, BENNU, D2R, Params, close_mission, is_bennu, make_cfg, optimize_transfer
from bennu_mission.core import TWO_PI, conic_points, time_of_flight
from bennu_mission.scenario2 import evaluate_transfer

REPO_URL = "https://github.com/sebastiano-golinelli/earth-to-bennu-mission-analysis"

# palette (kept in sync with .streamlit/config.toml)
BG, PANEL, GRID = "#060a13", "#0b1222", "rgba(150, 170, 210, 0.08)"
TEXT, MUTED = "#e7ecf6", "#8a97b1"
EARTH, TARGET, TRANSFER, SUN, VIOLET = "#4ea8ff", "#c9b48f", "#ff8a4c", "#ffd36b", "#9b7bff"
GOOD, BAD = "#3ddc97", "#ff6b6b"
COST_SCALE = [[0, "#fff1a8"], [0.12, "#ffc35a"], [0.3, "#ff8a4c"], [0.55, "#d6406b"],
              [0.8, "#6a2384"], [1, "#1b1035"]]                     # bright = cheap
ANGLE_SCALE = [[0, EARTH], [0.5, VIOLET], [1, TRANSFER]]
PLOT_CONFIG = {"scrollZoom": False, "displaylogo": False,
               "modeBarButtonsToRemove": ["select2d", "lasso2d", "autoScale2d"]}
CAMERAS = {"3D": dict(eye=dict(x=1.2, y=-1.95, z=1.0)),
           "Top": dict(eye=dict(x=0, y=0, z=2.4), up=dict(x=0, y=1, z=0)),
           "Edge-on": dict(eye=dict(x=0, y=-2.5, z=0.02))}

@st.cache_data(show_spinner="Searching 72,000 candidate transfers...")
def optimize(p: Params):
    return optimize_transfer(p)


@st.cache_data(show_spinner="Closing the mission: launch orbit, escape and capture...")
def full_mission(p: Params):
    opt = optimize(p)
    return None if opt is None else close_mission(p, opt["best"])


# ---- styling helpers -------------------------------------------------------------------------

CSS = f"""
<style>
[data-testid="stAppViewContainer"] {{
  background:
    radial-gradient(1px 1px at 12% 18%, rgba(255,255,255,.55), transparent 60%),
    radial-gradient(1px 1px at 78% 12%, rgba(255,255,255,.45), transparent 60%),
    radial-gradient(1.5px 1.5px at 42% 66%, rgba(255,255,255,.35), transparent 60%),
    radial-gradient(1px 1px at 88% 74%, rgba(255,255,255,.4), transparent 60%),
    radial-gradient(1px 1px at 24% 88%, rgba(255,255,255,.3), transparent 60%),
    radial-gradient(1200px 600px at 85% -10%, rgba(255,138,76,.10), transparent 60%),
    radial-gradient(900px 500px at -10% 30%, rgba(78,168,255,.08), transparent 60%),
    {BG};
}}
[data-testid="stHeader"] {{ background: transparent; }}
.block-container {{ padding-top: 2.2rem; max-width: 1280px; }}
.hero-eyebrow {{ font-family: 'JetBrains Mono', monospace; font-size: .75rem; letter-spacing: .18em;
  text-transform: uppercase; color: {TRANSFER}; margin-bottom: .35rem; }}
.hero-title {{ font-family: 'Space Grotesk', sans-serif; font-weight: 700; font-size: clamp(2rem, 4.5vw, 3.1rem);
  line-height: 1.05; margin: 0; padding: 0; color: {TEXT}; }}
.hero-title span {{ background: linear-gradient(90deg, {EARTH}, {VIOLET} 45%, {TRANSFER});
  -webkit-background-clip: text; background-clip: text; color: transparent; }}
.hero-sub {{ color: {MUTED}; font-size: 1.02rem; max-width: 760px; margin: .7rem 0 1rem; line-height: 1.55; }}
.pills {{ display: flex; flex-wrap: wrap; gap: .45rem; margin-bottom: 1.4rem; }}
.pill {{ font-size: .78rem; color: #c7d2e8; border: 1px solid #22304d; background: rgba(13,20,36,.7);
  padding: .22rem .65rem; border-radius: 999px; }}
.pill b {{ color: {GOOD}; font-weight: 600; }}
.kpis {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(165px, 1fr)); gap: .75rem; margin: .4rem 0 1.2rem; }}
.kpi {{ background: linear-gradient(180deg, rgba(17,26,46,.92), rgba(11,18,34,.92)); border: 1px solid #1c2740;
  border-radius: .8rem; padding: .85rem 1rem .8rem; position: relative; overflow: hidden; }}
.kpi::before {{ content: ""; position: absolute; inset: 0 0 auto 0; height: 2px;
  background: var(--accent, {TRANSFER}); opacity: .85; }}
.kpi .label {{ font-size: .7rem; letter-spacing: .12em; text-transform: uppercase; color: {MUTED}; }}
.kpi .value {{ font-family: 'JetBrains Mono', monospace; font-size: 1.55rem; font-weight: 500; color: {TEXT};
  margin-top: .25rem; white-space: nowrap; }}
.kpi .value small {{ font-size: .85rem; color: {MUTED}; margin-left: .25rem; }}
.kpi .note {{ font-size: .78rem; color: {MUTED}; margin-top: .15rem; }}
.kpi .delta-good {{ color: {GOOD}; }} .kpi .delta-bad {{ color: {BAD}; }}
.section-note {{ color: {MUTED}; font-size: .92rem; margin: -.2rem 0 .8rem; max-width: 900px; }}
.timeline {{ display: flex; flex-wrap: wrap; align-items: stretch; gap: .5rem; margin: .6rem 0 1rem; }}
.phase {{ flex: 1 1 150px; background: rgba(13,20,36,.85); border: 1px solid #1c2740; border-radius: .8rem;
  padding: .75rem .9rem; border-top: 2px solid var(--accent); }}
.phase .name {{ font-family: 'Space Grotesk', sans-serif; font-weight: 600; color: {TEXT}; font-size: .95rem; }}
.phase .dv {{ font-family: 'JetBrains Mono', monospace; color: var(--accent); font-size: 1.15rem; margin-top: .2rem; }}
.phase .meta {{ color: {MUTED}; font-size: .78rem; margin-top: .15rem; line-height: 1.4; }}
.arrow {{ align-self: center; color: #3a4a6b; font-size: 1.2rem; }}
.footer {{ color: {MUTED}; font-size: .82rem; border-top: 1px solid #1c2740; margin-top: 2.5rem; padding-top: 1rem; }}
.footer a, .pill a {{ color: #7cc0ff; }}
[data-testid="stTabs"] button p {{ font-family: 'Space Grotesk', sans-serif; font-size: 1rem; }}
@media (max-width: 640px) {{ .arrow {{ display: none; }} .kpi .value {{ font-size: 1.3rem; }} }}
</style>
"""


def html(markup):
    st.markdown(markup, unsafe_allow_html=True)


def kpi(label, value, unit="", note="", accent=TRANSFER, note_class=""):
    unit = f"<small>{unit}</small>" if unit else ""
    note = f'<div class="note {note_class}">{note}</div>' if note else ""
    return (f'<div class="kpi" style="--accent:{accent}"><div class="label">{label}</div>'
            f'<div class="value">{value}{unit}</div>{note}</div>')


def kpi_row(cards):
    html('<div class="kpis">' + "".join(cards) + "</div>")


def delta_note(value, ref, unit, scale=1.0):
    """Difference from the optimum as displayed (rounded), so the note never contradicts the numbers."""
    d = round(value * scale) - round(ref * scale)
    if d == 0:
        return "same as the optimum", "delta-good"
    return f"{d:+d} {unit} vs the optimum", "delta-bad" if d > 0 else "delta-good"


def styled(fig, height, title=None):
    fig.update_layout(
        height=height, margin=dict(l=64, r=16, t=52 if title else 12, b=56), showlegend=False,
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor=PANEL,
        font=dict(family="Inter, sans-serif", color="#c8d1e3", size=12),
        title=dict(text=title, font=dict(family="Space Grotesk, sans-serif", size=16, color=TEXT), x=0.01)
        if title else None,
        hoverlabel=dict(bgcolor="#121a2e", bordercolor="#2a3a5c", font=dict(family="JetBrains Mono, monospace", size=12)))
    fig.update_xaxes(gridcolor=GRID, zeroline=False, linecolor="#22304d", ticks="outside", tickcolor="#22304d",
                     automargin=True, title_standoff=8)
    fig.update_yaxes(gridcolor=GRID, zeroline=False, linecolor="#22304d", ticks="outside", tickcolor="#22304d",
                     automargin=True, title_standoff=8)
    return fig


def chart(fig):
    st.plotly_chart(fig, width="stretch", config=PLOT_CONFIG, theme=None)


def star(x, y, label):
    return go.Scatter(x=[x], y=[y], mode="markers", hovertemplate=f"{label}<extra></extra>",
                      marker=dict(symbol="star", size=17, color="white", line=dict(color=BG, width=1.5)))


# ---- plots ------------------------------------------------------------------------------------

def orbit_figure(cfg, sol, z_scale=1.0, view="3D"):
    nu = np.linspace(0, TWO_PI, 600)
    oE, oT = cfg.sc2.earth, cfg.sc2.target.orbit
    k = np.array([1, 1, z_scale])
    PE = conic_points(oE.a, oE.e, oE.i, oE.raan, oE.argp, nu) / AU * k
    PT = conic_points(oT.a, oT.e, oT.i, oT.raan, oT.argp, nu) / AU * k
    fig = go.Figure()
    for P, name, color in ((PE, "Earth orbit", EARTH), (PT, f"{cfg.sc2.target.name} orbit", TARGET)):
        fig.add_scatter3d(x=P[:, 0], y=P[:, 1], z=P[:, 2], mode="lines", name=name, hoverinfo="name",
                          line=dict(color=color, width=3))
    if sol.valid:
        nu_arc = sol.nu1_t + np.linspace(0, 1, 400) * np.deg2rad(sol.arc_deg)
        PX = conic_points(sol.a_t, sol.e_t, sol.i_t, sol.raan_t, sol.argp_t, nu_arc) / AU * k
        days = time_of_flight(sol.a_t, sol.e_t, sol.nu1_t, nu_arc, cfg.const.mu_sun) / 86400
        days[-1] = sol.tof_days                                  # a full arc wraps to 0 in time_of_flight
        fig.add_scatter3d(x=PX[:, 0], y=PX[:, 1], z=PX[:, 2], mode="lines", name="Transfer",
                          line=dict(color=days, colorscale=[[0, TRANSFER], [1, SUN]], width=7),
                          customdata=days, hovertemplate="transfer, day %{customdata:.0f}<extra></extra>",
                          showlegend=False)
        fig.add_scatter3d(x=[None], y=[None], z=[None], mode="lines", name="Transfer",   # legend entry only
                          line=dict(color=TRANSFER, width=7))
    for size, alpha in ((30, 0.07), (18, 0.16), (9, 1.0)):       # the Sun, with a soft glow
        fig.add_scatter3d(x=[0], y=[0], z=[0], mode="markers", showlegend=False, hovertext="Sun",
                          hoverinfo="text" if alpha == 1 else "skip", marker=dict(size=size, color=SUN, opacity=alpha))
    for r, name, color in ((sol.r1, "Departure", EARTH), (sol.r2, "Arrival", TARGET)):
        p = r / AU * k
        fig.add_scatter3d(x=[p[0]], y=[p[1]], z=[p[2]], mode="markers+text", name=name, text=[name],
                          textposition="top center", textfont=dict(color=color, size=12), hoverinfo="name",
                          showlegend=False, marker=dict(size=6, color=color, line=dict(color="white", width=1)))
    axis = dict(showbackground=False, gridcolor="rgba(150,170,210,0.07)", zerolinecolor="rgba(150,170,210,0.14)",
                color=MUTED, tickfont=dict(size=10))
    fig.update_layout(
        height=580, margin=dict(l=0, r=0, t=0, b=30), paper_bgcolor="rgba(0,0,0,0)",
        font=dict(family="Inter, sans-serif", color="#c8d1e3"),
        legend=dict(orientation="h", y=0.0, yanchor="top", x=0.5, xanchor="center", bgcolor="rgba(0,0,0,0)",
                    font=dict(size=12)),
        hoverlabel=dict(bgcolor="#121a2e", font=dict(family="JetBrains Mono, monospace")),
        scene=dict(aspectmode="data", camera=CAMERAS[view], xaxis=dict(axis, title="x [AU]"),
                   yaxis=dict(axis, title="y [AU]"),
                   zaxis=dict(axis, title="z [AU]" if z_scale == 1 else f"z [AU] × {z_scale:g}")))
    return fig


def cost_map_figure(opt, penalty):
    best = opt["best"]
    zmin = np.nanmin(opt["cost"])
    # angles are periodic: repeat the first row and column at 360 deg so the map closes
    cost = np.vstack([opt["cost"], opt["cost"][:1]])
    cost = np.hstack([cost, cost[:, :1]])
    fig = go.Figure(go.Contour(
        x=np.append(opt["vecs"][1], 360), y=np.append(opt["vecs"][0], 360), z=cost, colorscale=COST_SCALE, zmin=zmin, zmax=zmin + 8,
        contours=dict(start=zmin, end=zmin + 8, size=0.5, coloring="fill", showlines=True),
        line=dict(width=0.4, color="rgba(6,10,19,0.55)"), connectgaps=False,
        colorbar=dict(title=dict(text="cost [km/s]" if penalty else "Δv [km/s]", side="right"), thickness=12,
                      outlinewidth=0, tickmode="array", tickvals=np.arange(np.ceil(zmin), zmin + 8.01),
                      tickformat=".0f"),
        hovertemplate="arrival %{x:.0f}°<br>departure %{y:.0f}°<br>%{z:.3f} km/s<extra></extra>"))
    fig.add_trace(star(np.rad2deg(best.x[1]) % 360, np.rad2deg(best.x[0]) % 360, "optimum"))
    fig.update_xaxes(title="Arrival point on the target orbit [deg]", range=[0, 360], dtick=60)
    fig.update_yaxes(title="Departure point on Earth's orbit [deg]", range=[0, 360], dtick=60)
    return styled(fig, 440, "Cost map (Δv + time penalty)" if penalty else "Δv map: where to leave, where to arrive")


def tradeoff_figure(opt):
    best = opt["best"]
    keep = opt["dv"] <= best.dv + 5
    tof, dv, arc = opt["tof"][keep], opt["dv"][keep], opt["arc"][keep]
    if tof.size > 15000:                                       # keep the browser responsive
        idx = np.random.default_rng(0).choice(tof.size, 15000, replace=False)
        tof, dv, arc = tof[idx], dv[idx], arc[idx]
    fig = go.Figure(go.Scattergl(
        x=tof, y=dv, mode="markers", hovertemplate="TOF %{x:.0f} d<br>Δv %{y:.3f} km/s<extra></extra>",
        marker=dict(size=4, color=arc, colorscale=ANGLE_SCALE, opacity=0.75,
                    colorbar=dict(title=dict(text="transfer<br>angle [°]"), thickness=12, outlinewidth=0))))
    fig.add_trace(star(best.tof_days, best.dv, "optimum"))
    fig.update_xaxes(title="Time of flight [days]", rangemode="tozero")
    fig.update_yaxes(title="Δv departure + arrival [km/s]")
    return styled(fig, 440, "Cost versus time of flight")


def launch_window_figure(f):
    fig = go.Figure(go.Scatter(x=f["sweep_raan"], y=f["sweep_dv"], mode="lines", line=dict(color=EARTH, width=2.5),
                               fill="tozeroy", fillcolor="rgba(78,168,255,0.08)",
                               hovertemplate="RAAN %{x:.0f}°<br>Δv %{y:.3f} km/s<extra></extra>"))
    fig.add_trace(star(f["launch_raan"], f["sc1_dv"], "selected"))
    fig.update_xaxes(title="RAAN of the launch orbit [deg], set by the launch time", range=[0, 360], dtick=60)
    fig.update_yaxes(title="GTO → parking orbit Δv [km/s]", rangemode="tozero")
    return styled(fig, 360, "Launch window: first-phase cost over a day of launch times")


def budget_figure(f):
    vals = (f["sc1_dv"], f["escape_3d"], f["capture"], f["total"])
    fig = go.Figure(go.Waterfall(
        x=["GTO → parking", "Earth escape", "Capture", "Total"], y=list(vals), measure=["relative"] * 3 + ["total"],
        text=[f"{v:.3f}" for v in vals], textposition="outside",
        textfont=dict(family="JetBrains Mono, monospace", color=TEXT),
        connector=dict(line=dict(color="#2a3a5c", dash="dot")),
        increasing=dict(marker=dict(color=TRANSFER)), totals=dict(marker=dict(color=VIOLET)),
        hovertemplate="%{x}: %{y:.4f} km/s<extra></extra>"))
    fig.update_yaxes(title="Δv [km/s]", range=[0, f["total"] * 1.18])
    return styled(fig, 360, "Δv budget")


def why_invalid(sol, cfg):
    if not np.isfinite(sol.e_t):
        return "Departure and arrival points are aligned with the Sun: the transfer plane is undefined."
    if sol.e_t >= 1:
        return "Joining these two points with this orientation needs a hyperbola, not an ellipse."
    if sol.e_t < 0:
        return "No ellipse with this orientation passes through both points."
    if np.isfinite(sol.rp_t) and sol.rp_t <= cfg.sc2.min_perihelion:
        return (f"The transfer would pass {sol.rp_t / AU:.2f} AU from the Sun, closer than the "
                f"{cfg.sc2.min_perihelion / AU:.2f} AU limit.")
    return "This combination does not give a valid transfer."


# ---- page -------------------------------------------------------------------------------------

st.set_page_config(page_title="Earth → Bennu Mission Explorer", page_icon="🛰️", layout="wide")
html(CSS)

with st.sidebar:
    st.markdown("### Mission settings")
    st.caption("Every result on the page is recomputed when a setting changes.")
    with st.expander("🎯 Target orbit"):
        st.caption("Default: 101955 Bennu (JPL SBDB). A custom orbit keeps Bennu's mass and size.")
        if st.button("Reset to Bennu", width="stretch"):
            for key in ("a", "e", "i", "raan", "argp"):
                st.session_state.pop(f"tgt_{key}", None)
        a_au = st.number_input("Semi-major axis [AU]", 0.6, 3.0, BENNU.a / AU, 0.01, format="%.4f", key="tgt_a")
        ecc = st.number_input("Eccentricity", 0.0, 0.9, BENNU.e, 0.01, format="%.4f", key="tgt_e")
        inc = st.number_input("Inclination [deg]", 0.0, 60.0, BENNU.i / D2R, 0.5, format="%.3f", key="tgt_i")
        raan = st.number_input("RAAN [deg]", 0.0, 360.0, BENNU.raan / D2R, 1.0, format="%.3f", key="tgt_raan")
        argp = st.number_input("Argument of perihelion [deg]", 0.0, 360.0, BENNU.argp / D2R, 1.0, format="%.3f",
                               key="tgt_argp")
    st.markdown("**Constraints**")
    min_peri = st.slider("Closest approach to the Sun [AU]", 0.1, 0.9, 0.3, 0.05,
                         help="Transfers passing closer to the Sun are rejected (thermal limit).")
    long_arc = st.toggle("Allow transfers longer than half an orbit", True)
    tof_penalty = st.slider("Price of a day of flight [m/s]", 0.0, 20.0, 0.0, 0.5,
                            help="0 gives the cheapest transfer. A higher value trades Δv for a shorter trip.")
    st.markdown("**3D view**")
    view = st.segmented_control("Camera", list(CAMERAS), default="3D", label_visibility="collapsed") or "3D"
    z_scale = st.select_slider("Out-of-plane exaggeration", [1, 5, 10, 20], 1,
                               help="Bennu's orbit is tilted by only 6°: stretch the vertical axis to see it.")

params = Params(a_au, ecc, inc, raan, argp, min_peri, long_arc, tof_penalty)
cfg = make_cfg(params)
target = cfg.sc2.target.name

html(f"""
<div class="hero-eyebrow">Preliminary mission design · patched conics</div>
<div class="hero-title">Earth → {'Bennu' if is_bennu(params) else 'custom target'} <span>Mission Explorer</span></div>
<div class="hero-sub">From a geostationary transfer orbit to the near-Earth asteroid 101955 Bennu, the target of
NASA's OSIRIS-REx: find the best transfer, design your own, and close the whole mission budget.</div>
<div class="pills"><span class="pill"><b>✓</b> MATLAB ↔ Python agreement 1e-14</span>
<span class="pill"><b>✓</b> 37 automated tests</span><span class="pill">JPL Small-Body Database orbit</span>
<span class="pill">Open source · <a href="{REPO_URL}" target="_blank">GitHub</a></span></div>
""")

opt = optimize(params)
if opt is None:
    st.error("No valid transfer exists with these settings. Lower the Sun-distance limit or change the target orbit.")
    st.stop()
best = opt["best"]
mission, mission_error = None, ""
try:
    mission = full_mission(params)
except (RuntimeError, ValueError) as err:
    mission_error = str(err)

if mission:
    kpi_row([
        kpi("Total mission Δv", f"{mission['total']:.3f}", "km/s", "GTO to orbit at the target"),
        kpi("Mission duration", f"{mission['days']:.0f}", "days", f"{mission['tof']:.0f} days of cruise", EARTH),
        kpi("Arrival speed", f"{mission['v_inf_arr']:.2f}", "km/s", "relative to the target", TARGET),
        kpi("Launch orbit RAAN", f"{mission['launch_raan']:.0f}", "°", "best time of day to launch", VIOLET),
    ])

tab_traj, tab_design, tab_budget, tab_method = st.tabs(
    ["Optimal trajectory", "Design your own", "Mission budget", "Method & validation"])

with tab_traj:
    st.markdown(f"#### Cheapest direct transfer to {target}")
    html('<div class="section-note">Drag to rotate, hover for values. The colour of the transfer arc shows the '
         "days since departure. Use the sidebar to switch to a top or edge-on view.</div>")
    kpi_row([
        kpi("Δv departure", f"{best.dv1:.3f}", "km/s", "leave Earth's orbit", EARTH),
        kpi("Δv arrival", f"{best.dv2:.3f}", "km/s", "match the target's orbit", TARGET),
        kpi("Δv transfer", f"{best.dv:.3f}", "km/s", "departure + arrival"),
        kpi("Time of flight", f"{best.tof_days:.0f}", "days", f"transfer angle {best.arc_deg:.0f}°", SUN),
    ])
    chart(orbit_figure(cfg, best, z_scale, view))
    left, right = st.columns(2, gap="medium")
    with left:
        chart(cost_map_figure(opt, tof_penalty))
    with right:
        chart(tradeoff_figure(opt))
    st.caption(f"{opt['n_valid']:,} of {opt['n_total']:,} grid points give a valid elliptic transfer. "
               "The 12 best distinct candidates are refined with a gradient-based optimizer (SLSQP); "
               "the white star is the optimum.")

with tab_design:
    st.markdown("#### Choose the transfer yourself")
    html('<div class="section-note">A two-impulse transfer is fixed by where the spacecraft leaves Earth\'s orbit, '
         "where it meets the target, and how the transfer ellipse is oriented. Positions are true anomalies, "
         "measured from each orbit's perihelion.</div>")

    def load_optimum():
        x = np.rad2deg(np.mod(best.x, TWO_PI))
        st.session_state.update(nu1=float(round(x[0]) % 360), nu2=float(round(x[1]) % 360),
                                w=float(round(x[2]) % 360))

    if "nu1" not in st.session_state:
        load_optimum()
    c1, c2, c3 = st.columns(3)
    nu1 = c1.slider("Departure point on Earth's orbit", 0.0, 359.0, step=1.0, format="%.0f°", key="nu1")
    nu2 = c2.slider("Arrival point on the target orbit", 0.0, 359.0, step=1.0, format="%.0f°", key="nu2")
    w = c3.slider("Orientation of the transfer ellipse", 0.0, 359.0, step=1.0, format="%.0f°", key="w")
    st.button("↺ Back to the optimum", on_click=load_optimum)
    sol = evaluate_transfer(np.deg2rad([nu1, nu2, w]), cfg, "Manual")
    if sol.valid:
        dv_note, dv_cls = delta_note(sol.dv, best.dv, "m/s", 1000)
        tof_note, tof_cls = delta_note(sol.tof_days, best.tof_days, "d")
        kpi_row([
            kpi("Δv transfer", f"{sol.dv:.3f}", "km/s", dv_note, TRANSFER, dv_cls),
            kpi("Δv departure", f"{sol.dv1:.3f}", "km/s", "leave Earth's orbit", EARTH),
            kpi("Δv arrival", f"{sol.dv2:.3f}", "km/s", "match the target's orbit", TARGET),
            kpi("Time of flight", f"{sol.tof_days:.0f}", "days", tof_note, SUN, tof_cls),
        ])
    else:
        st.warning(why_invalid(sol, cfg), icon="⚠️")
    chart(orbit_figure(cfg, sol, z_scale, view))

with tab_budget:
    st.markdown("#### From the launch orbit to orbit around the asteroid")
    if mission is None:
        st.error(f"The mission cannot be closed with these settings: {mission_error}")
    else:
        f = mission
        html(f"""
<div class="timeline">
  <div class="phase" style="--accent:{EARTH}"><div class="name">GTO → parking orbit</div>
    <div class="dv">{f['sc1_dv']:.3f} km/s</div><div class="meta">strategy {f['sc1_name']}<br>{f['sc1_hours']:.1f} h</div></div>
  <div class="arrow">→</div>
  <div class="phase" style="--accent:{TRANSFER}"><div class="name">Earth escape</div>
    <div class="dv">{f['escape_3d']:.3f} km/s</div><div class="meta">v∞ {f['v_inf_dep']:.3f} km/s<br>{f['escape_days']:.1f} days to leave the SOI</div></div>
  <div class="arrow">→</div>
  <div class="phase" style="--accent:{SUN}"><div class="name">Heliocentric cruise</div>
    <div class="dv">coast</div><div class="meta">{f['tof']:.0f} days<br>transfer angle {best.arc_deg:.0f}°</div></div>
  <div class="arrow">→</div>
  <div class="phase" style="--accent:{TARGET}"><div class="name">Capture at the target</div>
    <div class="dv">{f['capture']:.3f} km/s</div><div class="meta">{f['capture_alt']:g} km orbit<br>SOI only {f['soi_target']:.2f} km</div></div>
</div>""")
        left, right = st.columns([1, 1.25], gap="medium")
        with left:
            chart(budget_figure(f))
        with right:
            chart(launch_window_figure(f))
        with st.expander("Geocentric strategies compared (S = shape, P = plane, W = periapsis; a/p = apsides used)"):
            st.dataframe(f["sc1_table"], hide_index=True, width="stretch")
        st.caption(f"The escape is computed in full 3D from the parking orbit; the coplanar estimate gives "
                   f"{f['escape_2d']:.4f} km/s, which confirms that the parking orbit is aligned with the departure "
                   "asymptote. At the target the 'capture' is in practice a rendezvous: the sphere of influence is "
                   "only a few kilometres wide.")

with tab_method:
    st.markdown(f"""
#### What is computed
A preliminary mission analysis with **patched conics** and **impulsive manoeuvres**, in three phases:
1. **Geocentric transfer** from a GTO (Cape Canaveral, 185 × 35 786 km, 27°) to a parking orbit designed so that a
   single tangential burn leaves Earth exactly along the required departure direction. 16 manoeuvre strategies are
   compared, and the launch time (orientation of the GTO) is optimized.
2. **Heliocentric transfer**: a two-impulse ellipse from Earth's orbit to the target's orbit, searched over
   72,000 combinations and refined with SLSQP.
3. **Escape and capture** with hyperbolic arcs inside the spheres of influence.

#### How it is validated
The same analysis exists in **MATLAB** and **Python**: the Python results match MATLAB to a relative error of 1e-14.
Both versions run automated tests on every commit (GitHub Actions), including numerical propagation of the transfer
and the escape asymptote recomputed from the post-burn state. This app reuses the tested package unchanged; its
optimizer was cross-checked against the full global search (differential evolution and 200-point multi-start).

#### Limits
Two-body arcs, impulsive burns, mean J2000 Earth orbit and no launch date: the result is the best *geometry*, not a
calendar launch window. No gravity assists. Good for a first sizing, not for flight.

**Stack:** Python · NumPy · SciPy · Plotly · Streamlit. Source code on [GitHub]({REPO_URL}).
""")

html(f"""<div class="footer">Built by Sebastiano Enrico Golinelli, aerospace engineering student ·
<a href="{REPO_URL}" target="_blank">source code</a> · orbit data: JPL Small-Body Database ·
Bennu GM: OSIRIS-REx radio science</div>""")
