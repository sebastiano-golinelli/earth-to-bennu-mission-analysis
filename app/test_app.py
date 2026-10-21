"""Tests of the web app and of the computations behind it.

Run from the repository root:  pytest app
"""

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from explorer import Params, close_mission, make_cfg, optimize_transfer
from bennu_mission.scenario2 import run_scenario2

APP = str(Path(__file__).with_name("streamlit_app.py"))


def page_text(at):
    return " ".join(m.value for m in at.markdown)


def test_page_shows_the_validated_bennu_results():
    at = AppTest.from_file(APP, default_timeout=180)
    at.run()
    assert not at.exception
    text = page_text(at)
    assert "5.087" in text            # heliocentric transfer [km/s]
    assert "5.660" in text            # total mission delta-v [km/s]
    assert ">123<" in text            # mission duration [days]


def test_infeasible_manual_transfer_is_explained():
    at = AppTest.from_file(APP, default_timeout=180)
    at.run()
    at.slider(key="nu1").set_value(10.0)
    at.slider(key="nu2").set_value(200.0)
    at.run()
    assert not at.exception
    assert at.warning, "an infeasible transfer must show why"


def test_mission_budget_matches_run_mission():
    # reference values of python/tests (run_mission with the full global search)
    p = Params()
    best = optimize_transfer(p)["best"]
    m = close_mission(p, best)
    assert best.dv == pytest.approx(5.086991961, abs=1e-6)
    assert m["total"] == pytest.approx(5.660054857, abs=1e-5)
    assert m["days"] == pytest.approx(122.73, abs=0.05)
    assert m["launch_raan"] == pytest.approx(192)


@pytest.mark.parametrize("orbit", [
    (2.1388, 0.1421, 16.2137, 250.6378, 105.3795),    # refining only the best grid point stops at 11.24 km/s
    (0.9199, 0.1184, 20.0319, 209.5783, 33.8863),     # ... and here at 10.57 km/s
    (1.4064, 0.2395, 3.9935, 264.4478, 40.9219),
])
def test_optimizer_matches_the_full_global_search(orbit):
    p = Params(*orbit)
    cfg = make_cfg(p)
    cfg.sc2.de.seeds = [1, 2, 3]
    reference = run_scenario2(cfg).best.J
    assert optimize_transfer(p)["best"].J == pytest.approx(reference, abs=1e-4)
