function res = runMission(cfg)
%RUNMISSION Run the three scenarios in sequence and print a summary.
%   res = RUNMISSION(cfg) returns the results of each scenario (res.sc1,
%   res.sc2, res.sc3) and the mission budget (res.mission). Plots are drawn
%   when cfg.plot.enable is true.

% The heliocentric transfer fixes the departure asymptote, which in turn
% defines the parking orbit that the geocentric transfer has to reach.
res.sc2 = runScenario2(cfg);
vInf = eclipticToEquatorial(res.sc2.best.vT1 - res.sc2.best.vEarth, cfg.const.obliquity);
res.sc1 = runScenario1(cfg, vInf);
res.sc3 = runScenario3(cfg, res.sc1, res.sc2);

m.dvGeocentric = res.sc1.best.dvTotal;
m.dvCoplanar = m.dvGeocentric + res.sc3.dvCoplanar;
m.dvNonCoplanar = m.dvGeocentric + res.sc3.dvNonCoplanar;
m.durationDays = res.sc1.best.time/86400 + res.sc3.escape.timeInSoi/86400 + ...
    res.sc2.best.tofDays + res.sc3.capture.timeInSoi/86400;
res.mission = m;

if cfg.verbose
    printSummary(cfg, res);
end
if cfg.plot.enable
    plotScenario1(cfg, res.sc1);
    plotScenario2(cfg, res.sc2);
    plotScenario3(cfg, res.sc3);
end
end

function printSummary(cfg, res)
r1 = res.sc1; r2 = res.sc2; r3 = res.sc3; m = res.mission;
line = repmat('=', 1, 72);

fprintf('\n%s\nSCENARIO 1 - Geocentric transfer: GTO -> parking orbit\n%s\n', line, line);
if ~isempty(r1.launchSweep)
    fprintf('Launch RAAN sweep: geocentric dv from %.3f to %.3f km/s -> launch RAAN = %.0f deg\n', ...
        min(r1.launchSweep.DeltaV_kms), max(r1.launchSweep.DeltaV_kms), rad2deg(r1.initial.raan));
end
fprintf('Parking orbit (%s): a = %.3f km, e = %.6f, i = %.4f deg, RAAN = %.4f deg, argp = %.4f deg\n', ...
    cfg.sc1.parking.mode, r1.final.a, r1.final.e, rad2deg([r1.final.i, r1.final.raan, r1.final.argp]));
if isfield(r1.final, 'asymptote')
    fprintf('  departure asymptote: RA = %.3f deg, DEC = %.3f deg\n', ...
        rad2deg(r1.final.asymptote.rightAscension), rad2deg(r1.final.asymptote.declination));
end
disp(r1.table);
fprintf('Selected (%s): %s, dv = %.4f km/s, time = %.2f h\n', cfg.sc1.selection, ...
    r1.best.name, r1.best.dvTotal, r1.best.time/3600);

fprintf('\n%s\nSCENARIO 2 - Heliocentric transfer: Earth -> %s\n%s\n', line, cfg.sc2.target.name, line);
fprintf('Valid grid points: %d / %d  (long arcs allowed: %d)\n', nnz(r2.gridRecords.valid), ...
    numel(r2.gridRecords.valid), cfg.sc2.allowLongArc);
disp(r2.comparison);
if ~isempty(r2.localMinima)
    fprintf('MultiStart: %d distinct valid local minima; best five:\n', height(r2.localMinima));
    disp(r2.localMinima(1:min(5, end), :));
end
b = r2.best;
fprintf('Best (%s): dv = %.4f km/s (%.4f + %.4f), TOF = %.1f days, arc = %.1f deg\n', ...
    b.method, b.dv, b.dv1, b.dv2, b.tofDays, b.arcDeg);
fprintf('  transfer orbit: a = %.4f AU, e = %.4f, i = %.3f deg\n', b.aT/cfg.const.AU, b.eT, rad2deg(b.iT));

fprintf('\n%s\nSCENARIO 3 - Patched conics\n%s\n', line, line);
fprintf('SOI radius: Earth %.0f km, %s %.3f km\n', r3.soiEarth, cfg.sc2.target.name, r3.soiTarget);
fprintf('Escape (coplanar, tangential at parking periapsis): v_inf = %.4f km/s, dv = %.4f km/s, e = %.4f\n', ...
    r3.vInfDep, r3.escape.dv, r3.escape.e);
fprintf('Escape (non-coplanar, best burn point nu = %.1f deg): dv = %.4f km/s\n', ...
    rad2deg(r3.escape3D.nuBurnOnParking), r3.escape3D.dv);
fprintf('Capture: v_inf = %.4f km/s\n', r3.vInfArrMag);
disp(r3.captureTable);

fprintf('\n%s\nMISSION BUDGET\n%s\n', line, line);
fprintf('Geocentric transfer (%s)          %8.4f km/s\n', r1.best.name, m.dvGeocentric);
fprintf('Escape + capture, coplanar escape      %8.4f km/s\n', r3.dvCoplanar);
fprintf('Escape + capture, non-coplanar escape  %8.4f km/s\n', r3.dvNonCoplanar);
fprintf('TOTAL (coplanar, idealized)            %8.4f km/s\n', m.dvCoplanar);
fprintf('TOTAL (non-coplanar, realistic)        %8.4f km/s\n', m.dvNonCoplanar);
fprintf('Mission duration                       %8.1f days\n', m.durationDays);
end
