function exportMatlabReference()
%EXPORTMATLABREFERENCE Save the MATLAB results used by the Python parity tests.
%   Runs the MATLAB mission and writes matlab_reference.json in this folder.
%   Run it again after any change to the MATLAB model.

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fullfile(root, 'config'));
addpath(genpath(fullfile(root, 'src')));
cfg = missionConfig();
cfg.plot.enable = false;
cfg.verbose = false;
r = runMission(cfg);

f = r.sc1.final;
ref.sc1 = struct('launchRaanDeg', rad2deg(r.sc1.initial.raan), ...
    'launchSweepDv', r.sc1.launchSweep.DeltaV_kms, ...
    'strategies', r.sc1.table.Strategy, 'dv', r.sc1.table.DeltaV_kms, 'timeH', r.sc1.table.Time_h, ...
    'selected', r.sc1.best.name, 'parking', [f.a f.e f.i f.raan f.argp]);

rec = r.sc2.gridRecords;
idx = (1:97:numel(rec.dv)).';
b = r.sc2.best;
ref.sc2 = struct('gridValidCount', nnz(rec.valid), 'gridSampleIndex', idx - 1, ...
    'gridSampleDv', rec.dv(idx), 'gridBestX', r.sc2.grid.x, 'gridBestDv', r.sc2.grid.dv, ...
    'fminconDv', r.sc2.fmincon.dv, 'gaRunsJ', r.sc2.gaRuns(:, 2), ...
    'multiStartDv', r.sc2.multiStart.dv, 'localMinima', height(r.sc2.localMinima), ...
    'bestX', b.x, 'bestDv', b.dv, 'bestDv1', b.dv1, 'bestDv2', b.dv2, ...
    'bestTofDays', b.tofDays, 'bestArcDeg', b.arcDeg);

s3 = r.sc3;
ref.sc3 = struct('soiEarth', s3.soiEarth, 'soiTarget', s3.soiTarget, ...
    'escapeDv', s3.escape.dv, 'escapeEcc', s3.escape.e, 'escapeTimeInSoi', s3.escape.timeInSoi, ...
    'escape3dDv', s3.escape3D.dv, 'burnNuDeg', rad2deg(s3.escape3D.nuBurnOnParking), ...
    'burnSweepDv', s3.burnSweep.dv(:), 'captureDvCircular', s3.captureTable.DeltaVCircular_kms, ...
    'captureSrpRatio', s3.captureTable.SRPRatio, 'dvCoplanar', s3.dvCoplanar, ...
    'dvNonCoplanar', s3.dvNonCoplanar);

m = r.mission;
ref.mission = struct('dvGeocentric', m.dvGeocentric, 'dvCoplanar', m.dvCoplanar, ...
    'dvNonCoplanar', m.dvNonCoplanar, 'durationDays', m.durationDays);
ref.meta = struct('matlabVersion', version, 'created', char(datetime('now', 'Format', 'yyyy-MM-dd')));

fid = fopen(fullfile(here, 'matlab_reference.json'), 'w');
fwrite(fid, jsonencode(ref, 'PrettyPrint', true));
fclose(fid);
fprintf('Saved %s\n', fullfile(here, 'matlab_reference.json'));
end
