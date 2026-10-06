import json
import sys
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bennu_mission import mission_config, run_mission  # noqa: E402
from bennu_mission.core import ecliptic_to_equatorial  # noqa: E402
from bennu_mission.scenario1 import run_scenario1  # noqa: E402
from bennu_mission.scenario2 import evaluate_transfer  # noqa: E402
from bennu_mission.scenario3 import run_scenario3  # noqa: E402


@pytest.fixture(scope="session")
def cfg():
    c = mission_config()
    c.plot.enable = False
    c.verbose = False
    return c


@pytest.fixture(scope="session")
def ref():
    """MATLAB results (reference/exportMatlabReference.m); null -> NaN."""
    def to_numpy(d):
        return {k: (to_numpy(v) if isinstance(v, dict)
                    else np.array([np.nan if x is None else x for x in v], float)
                    if isinstance(v, list) and v and not isinstance(v[0], str) else v)
                for k, v in d.items()}
    return to_numpy(json.loads((ROOT / "reference" / "matlab_reference.json").read_text()))


@pytest.fixture(scope="session")
def from_matlab_optimum(cfg, ref):
    """Deterministic chain started from the MATLAB optimal transfer: SC2 point -> SC1 -> SC3."""
    best = evaluate_transfer(ref["sc2"]["bestX"], cfg, "MATLAB optimum")
    r2 = NS(best=best)
    v_inf = ecliptic_to_equatorial(best.v_t1 - best.v_earth, cfg.const.obliquity)
    r1 = run_scenario1(cfg, v_inf)
    r3 = run_scenario3(cfg, r1, r2)
    return NS(sc1=r1, sc2=r2, sc3=r3)


@pytest.fixture(scope="session")
def full_run(cfg):
    """Complete Python run, optimizers included."""
    return run_mission(cfg)
