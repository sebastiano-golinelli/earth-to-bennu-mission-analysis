"""Earth -> Bennu mission analysis, Python version.

Run from the python/ folder:  python main.py
Edit bennu_mission/config.py to change data and options.
"""

import time

from bennu_mission import mission_config, run_mission

if __name__ == "__main__":
    cfg = mission_config()
    t0 = time.perf_counter()
    results = run_mission(cfg)
    print(f"\nRun time: {time.perf_counter() - t0:.1f} s")
