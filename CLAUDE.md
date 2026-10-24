# Earth → Bennu mission analysis

Preliminary mission analysis from a GTO to asteroid 101955 Bennu, in MATLAB with a full Python port and
a Streamlit web app.

## Layout
- `main.m`, `src/`, `config/`, `tests/` (15 `matlab.unittest` tests), `figures/`: MATLAB version
- `python/`: Python port (`bennu_mission/`, 16 pytest tests, `compare_with_matlab.py`, `reference/`)
- `app/`: Streamlit web app (`streamlit_app.py`, computations in `explorer.py`, 6 tests); theme in
  `.streamlit/config.toml` at the root
- `.github/workflows/`: CI with `matlab-tests.yml`, `python-tests.yml`, `app-tests.yml`

## Commands
- MATLAB: `runtests('tests')` (about 35 s)
- Python: `cd python`, activate `.venv`, then `pytest` and `python compare_with_matlab.py`
- App: `streamlit run app/streamlit_app.py`

## Reference results (must not change without a reason)
Total **5.660 km/s, ≈123 days**: GTO → SPW parking orbit 1.755 km/s in 49 h; 3D escape 0.683; arrival
3.223; 116-day cruise. Python and MATLAB agree to within 1e-13.

## Rules
- Every change must keep **all three CI workflows green**: run the tests before proposing a push.
- MATLAB and Python must keep matching the reference results above; `compare_with_matlab.py` checks it.
