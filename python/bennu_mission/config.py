"""All data and options of the Earth -> Bennu mission (same values as config/missionConfig.m).

Units: km, km/s, s and rad unless stated otherwise.
"""

from types import SimpleNamespace as NS

import numpy as np

from .core import Orbit

D2R = np.pi / 180


def mission_config():
    c = NS(
        mu_sun=1.32712440018e11,          # [km^3/s^2] Sun (JPL DE405)
        mu_earth=398600.4418,             # [km^3/s^2] Earth (WGS-84)
        r_earth=6378.137,                 # [km] Earth equatorial radius (WGS-84)
        AU=149597870.7,                   # [km] astronomical unit (IAU 2012)
        obliquity=84381.406 / 3600 * D2R,  # [rad] obliquity of the ecliptic at J2000 (IAU 2006)
    )

    # Scenario 1 - geocentric transfer: GTO from Cape Canaveral (185 x 35786 km, 27 deg)
    rp, ra = c.r_earth + 185, c.r_earth + 35786
    sc1 = NS(
        initial=Orbit(a=(rp + ra) / 2, e=(ra - rp) / (ra + rp), i=27 * D2R, raan=np.nan,
                      argp=178 * D2R, nu=20 * D2R),
        launch_raan="optimize",           # 'optimize' or a fixed value [deg]
        launch_raan_step=1,               # [deg]
        parking=NS(mode="designed",       # 'designed' (aligned with the asymptote) or 'given'
                   periapsis_altitude=600, apoapsis_altitude=60000,
                   state_r=np.array([28791.1999, -20951.2138, 22784.6716]),
                   state_v=np.array([2.679913, 0.068450, 0.896740])),
        strategies=["SPW", "PSW", "PWS", "S(P)W"],
        types=["ap", "pa", "pp", "aa"],
        selection="minDeltaV",            # 'minDeltaV', 'minTime' or a strategy name, e.g. 'SPW ap'
    )

    # Scenario 2 - heliocentric transfer Earth -> Bennu (ecliptic J2000)
    sc2 = NS(
        earth=Orbit(a=1.00000261 * c.AU, e=0.01671123, i=0.0, raan=0.0, argp=102.93768193 * D2R),
        target=NS(name="101955 Bennu", orbit=Orbit(a=1.126391025894812 * c.AU, e=0.2037450762416414,
                                                   i=6.03494377024794 * D2R, raan=2.06086619569642 * D2R,
                                                   argp=66.22306084084298 * D2R),
                  mu=4.8904e-9, radius=0.48444 / 2),
        allow_long_arc=True,
        min_perihelion=0.3 * c.AU,        # [km]
        weight_tof=0.0,                   # J = dv + weight_tof * TOF[days]
        grid=(30, 30, 80),
        local=NS(max_iterations=300),
        de=NS(seeds=[1, 2, 3, 4, 5], generations=80, population=81),
        multistart=NS(points=200, seed=7),
    )

    # Scenario 3 - patched conics
    sc3 = NS(
        soi_mode="semimajor",             # 'semimajor' or 'patch'
        capture_altitudes=np.array([0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2]),   # [km]
        nominal_capture_altitude=1.0,     # [km]
        escape_altitude_sweep=np.arange(300, 5001, 200),                       # [km]
        burn_sweep_points=361,
        srp=NS(Cr=1.5, area_to_mass=0.02, pressure_1au=4.56e-6),             # [-], [m^2/kg], [N/m^2]
    )

    plot = NS(enable=True, export=False, folder="figures")
    return NS(const=c, sc1=sc1, sc2=sc2, sc3=sc3, verbose=True, plot=plot)
