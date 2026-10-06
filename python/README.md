# Python version

A complete port of the MATLAB analysis to Python (NumPy, SciPy, Matplotlib), checked against the MATLAB results with an automated test suite.

## Results: Python versus MATLAB

Where the computation is deterministic, the two versions are compared starting from the same point (the MATLAB optimal transfer), so the comparison isolates the ported formulas. They agree to machine precision:

| Quantity | Values compared | Max relative difference |
|---|---|---|
| Launch window sweep (Δv) | 360 | 9.5 × 10⁻¹⁵ |
| 16 geocentric strategies (Δv and time) | 32 | 7.5 × 10⁻¹⁶ |
| Parking orbit elements | 5 | 6.9 × 10⁻¹⁶ |
| Heliocentric grid search (sample of Δv) | 321 | 8.9 × 10⁻¹⁴ |
| Optimal transfer (Δv, TOF) | 2 | 1.7 × 10⁻¹⁶ |
| 3D escape, burn point sweep (Δv) | 252 | 4.8 × 10⁻¹⁴ |
| Capture table (Δv) | 8 | 2.8 × 10⁻¹⁶ |

The optimizers are different algorithms, so only their final optimum is compared:

| Python | MATLAB | Δv difference [km/s] |
|---|---|---|
| SLSQP | `fmincon` (SQP) | 1.6 × 10⁻¹² |
| Differential evolution + SLSQP | `ga` + `fmincon` | 5.2 × 10⁻¹³ |
| Multi-start SLSQP | `MultiStart` | 5.1 × 10⁻¹⁴ |

The figures in [`figures/`](figures) are the Matplotlib counterparts of the MATLAB ones.

## How to run

```bash
cd python
python -m venv .venv
.venv\Scripts\activate          # Windows (on Linux/macOS: source .venv/bin/activate)
pip install -r requirements.txt

python main.py                  # full mission, summary and figures (about 30 s)
pytest                          # 16 tests, including the MATLAB comparison
python compare_with_matlab.py   # prints the comparison tables above
```

The MATLAB reference values are in `reference/matlab_reference.json`, created by `reference/exportMatlabReference.m`.

## Porting notes

| MATLAB | Python | Note |
|---|---|---|
| structs | `dataclass` (orbits, results), `SimpleNamespace` (configuration) | same dot notation as MATLAB |
| column vectors 3×1 | arrays of shape `(3,)`, stacks `(N, 3)` | components on the last axis |
| `for` loops over the grid | vectorized `evaluate_transfers` | 72,000 transfers in 0.07 s instead of 6.1 s |
| `fmincon` (SQP) | `scipy.optimize.minimize`, SLSQP | gradient by forward differences in one vectorized call |
| `ga` | `scipy.optimize.differential_evolution` | whole population evaluated at once (`vectorized=True`) |
| `MultiStart` | SLSQP from 200 random points | no direct SciPy equivalent |
| `fzero`, `ode113` | `brentq`, `solve_ivp` (DOP853) | |
| `matlab.unittest` | `pytest` | |

Lessons learned during the port:

- **Vectorizing is not enough.** The first version evaluated one transfer in about 1 ms because of the overhead of building small NumPy arrays: the whole run took 73 s. Writing the rotation matrices in closed form and computing the gradient in a single vectorized call brought it down to about 25–30 s, the same as MATLAB.
- **Same optimum, different paths.** SLSQP takes two trial steps into the invalid region of the search space (where the objective is a penalty) before converging, while `fmincon` does not. With the same budget, differential evolution reaches the global optimum with 1 seed out of 5, plus 1 seed in a neighbouring minimum only 0.6 m/s higher; `ga` reaches it with 1 seed out of 5.
- **Plot autoscaling differs.** With a fixed y range, MATLAB scales the x axis on the visible points only, Matplotlib on all of them: a few valid but extreme transfers (times of flight of millions of days, Δv above 25 km/s) stretched the axis until the limits were set explicitly.

## Structure

```
main.py                    entry point
compare_with_matlab.py     comparison tables
bennu_mission/
    config.py              data and options (same values as config/missionConfig.m)
    core.py                Keplerian elements, frames, time of flight
    maneuvers.py           bitangent transfer, plane change, periapsis rotation
    scenario1.py           geocentric strategies, parking orbit design, launch window
    scenario2.py           vectorized transfer model and optimizers
    scenario3.py           escape and capture with patched conics
    mission.py             runs the three phases in sequence
    plotting.py            figures
tests/                     pytest suite (core checks and MATLAB comparison)
reference/                 MATLAB reference results
```
