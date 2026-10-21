"""Smoke test of the web app: every tab renders and shows the validated results.

Run from the repository root:  pytest app
"""

from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).with_name("streamlit_app.py"))


def metrics(at):
    return {m.label: m.value for m in at.metric}


def test_optimal_transfer_and_full_mission():
    at = AppTest.from_file(APP, default_timeout=180)
    at.run()
    assert not at.exception
    assert metrics(at)["Total Δv"] == "5.087 km/s"

    next(b for b in at.button if "full mission" in b.label).click()
    at.run()
    assert not at.exception
    m = metrics(at)
    assert m["Total Δv"] == "5.660 km/s"
    assert m["Mission duration"] == "123 days"


def test_infeasible_manual_transfer_is_explained():
    at = AppTest.from_file(APP, default_timeout=180)
    at.run()
    at.slider(key="nu1").set_value(10.0)
    at.slider(key="nu2").set_value(200.0)
    at.run()
    assert not at.exception
    assert at.warning, "an infeasible transfer must show why"
