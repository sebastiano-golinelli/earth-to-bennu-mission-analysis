function cfg = missionConfig()
%MISSIONCONFIG All data and options of the Earth -> Bennu mission analysis.
%   Units: km, km/s, s and rad unless stated otherwise.

d2r = pi/180;

%% Constants
c.muSun   = 1.32712440018e11;        % [km^3/s^2] Sun (JPL DE405)
c.muEarth = 398600.4418;             % [km^3/s^2] Earth (WGS-84)
c.rEarth  = 6378.137;                % [km] Earth equatorial radius (WGS-84)
c.AU      = 149597870.7;             % [km] astronomical unit (IAU 2012)
c.obliquity = 84381.406/3600*d2r;    % [rad] obliquity of the ecliptic at J2000 (IAU 2006)
cfg.const = c;

%% Scenario 1 - geocentric transfer: launch orbit -> parking orbit
% Launch orbit: GTO from Cape Canaveral (185 x 35786 km, 27 deg), shortly
% after injection.
rp = c.rEarth + 185;
ra = c.rEarth + 35786;
cfg.sc1.initial = struct('a', (rp + ra)/2, 'e', (ra - rp)/(ra + rp), 'i', 27*d2r, ...
    'raan', NaN, 'argp', 178*d2r, 'nu', 20*d2r);
% The RAAN of the launch orbit is set by the launch time: 'optimize' sweeps
% it and keeps the cheapest geocentric transfer, a number fixes it [deg].
cfg.sc1.launchRaan = 'optimize';
cfg.sc1.launchRaanStep = 1;          % [deg]
% Parking orbit from which the escape burn is performed.
%  'designed': highly elliptical orbit whose plane contains the departure
%              asymptote of Scenario 2 and whose periapsis is placed for a
%              tangential escape burn (see designParkingOrbit).
%  'given'   : orbit given as a state vector (arrival point = that state).
cfg.sc1.parking.mode = 'designed';
cfg.sc1.parking.periapsisAltitude = 600;      % [km]
cfg.sc1.parking.apoapsisAltitude = 60000;     % [km]
cfg.sc1.parking.state.r = [28791.1999; -20951.2138; 22784.6716];   % [km]   ('given' mode)
cfg.sc1.parking.state.v = [2.679913; 0.068450; 0.896740];          % [km/s] ('given' mode)
cfg.sc1.strategies = {'SPW', 'PSW', 'PWS', 'S(P)W'};
cfg.sc1.types = {'ap', 'pa', 'pp', 'aa'};
cfg.sc1.selection = 'minDeltaV';     % 'minDeltaV', 'minTime' or a strategy name, e.g. 'SPW ap'

%% Scenario 2 - heliocentric transfer: Earth -> Bennu
% Earth: mean J2000 elements (Standish, JPL "Approximate Positions of the
% Planets"); ecliptic J2000 frame, inclination set to 0 by definition.
cfg.sc2.earth = struct('name', "Earth", 'a', 1.00000261*c.AU, 'e', 0.01671123, ...
    'i', 0, 'raan', 0, 'argp', 102.93768193*d2r);
% 101955 Bennu: JPL Small-Body Database, orbit solution 118 (epoch
% 2011-01-01 TDB, ecliptic J2000). GM from Chesley et al. 2020, diameter
% from Daly et al. 2020 (OSIRIS-REx data).
cfg.sc2.target = struct('name', "101955 Bennu", 'a', 1.126391025894812*c.AU, ...
    'e', 0.2037450762416414, 'i', 6.03494377024794*d2r, 'raan', 2.06086619569642*d2r, ...
    'argp', 66.22306084084298*d2r, 'mu', 4.8904e-9, 'radius', 0.48444/2);
% Transfer model
cfg.sc2.allowLongArc = true;         % false: only arcs < 180 deg (original lab model)
cfg.sc2.minPerihelion = 0.3*c.AU;    % [km] thermal limit on the transfer perihelion
cfg.sc2.weightTOF = 0;               % J = dv + weightTOF*TOF[days]
% Optimizers
cfg.sc2.grid = [30 30 80];           % points in nuDep, nuArr, argpT
cfg.sc2.fmincon.maxIterations = 300;
cfg.sc2.ga.seeds = 1:5;
cfg.sc2.ga.generations = 80;
cfg.sc2.ga.population = 80;
cfg.sc2.multiStart.points = 200;     % 0 disables MultiStart
cfg.sc2.multiStart.seed = 7;

%% Scenario 3 - patched conics: Earth escape and capture at Bennu
cfg.sc3.soiMode = 'semimajor';       % 'semimajor' (classical) or 'patch' (distance at the patch point)
% Bennu's sphere of influence is only ~2.8 km: capture radii must stay
% well inside it. OSIRIS-REx orbited at radii of roughly 1-2 km.
cfg.sc3.captureAltitudes = [0.05 0.1 0.25 0.5 0.75 1 1.5 2];   % [km] above the mean radius
cfg.sc3.nominalCaptureAltitude = 1;                            % [km]
cfg.sc3.escapeAltitudeSweep = 300:200:5000;                    % [km] parking periapsis altitudes
cfg.sc3.burnSweepPoints = 361;                                 % burn points along the parking orbit
% Spacecraft model for the solar radiation pressure check
cfg.sc3.srp = struct('Cr', 1.5, 'areaToMass', 0.02, 'pressure1AU', 4.56e-6);   % [-], [m^2/kg], [N/m^2]

%% Output
cfg.verbose = true;                  % print the summary tables
cfg.plot.enable = true;
cfg.plot.allSC1Cases = false;        % true: one figure for each of the 16 strategies
cfg.plot.export = false;             % true: save PNG files in cfg.plot.folder
cfg.plot.folder = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
end
