function r1 = runScenario1(cfg, vInf)
%RUNSCENARIO1 Geocentric transfer from the launch orbit to the parking orbit.
%   r1 = RUNSCENARIO1(cfg, vInf) evaluates every combination of strategy
%   (SPW, PSW, PWS, S(P)W) and bitangent type (ap, pa, pp, aa) from the
%   launch orbit cfg.sc1.initial to the parking orbit, and selects one
%   according to cfg.sc1.selection.
%   The parking orbit is designed around the departure excess velocity vInf
%   (3x1, equatorial frame) or given as a state vector (cfg.sc1.parking).
%   With cfg.sc1.launchRaan = 'optimize', the RAAN of the launch orbit (set
%   by the launch time) is swept and the cheapest one is used.

mu = cfg.const.muEarth;
oI = cfg.sc1.initial;

if ischar(cfg.sc1.launchRaan) && strcmp(cfg.sc1.launchRaan, 'optimize')
    raan = (0:cfg.sc1.launchRaanStep:360 - cfg.sc1.launchRaanStep).';
    dv = nan(size(raan)); time = dv; name = strings(size(raan));
    for k = 1:numel(raan)
        oI.raan = deg2rad(raan(k));
        [cases, T] = evaluateAll(cfg, oI, parkingOrbit(cfg, vInf, oI), mu);
        kSel = selectCase(T, cfg.sc1.selection);
        dv(k) = cases(kSel).dvTotal; time(k) = T.Time_h(kSel); name(k) = cases(kSel).name;
    end
    r1.launchSweep = table(raan, dv, time, name, ...
        'VariableNames', {'LaunchRAAN_deg', 'DeltaV_kms', 'Time_h', 'Strategy'});
    [~, kBest] = min(dv);
    oI.raan = deg2rad(raan(kBest));
else
    oI.raan = deg2rad(cfg.sc1.launchRaan);
    r1.launchSweep = table();
end

[r1.initialState.r, r1.initialState.v] = kep2car(oI.a, oI.e, oI.i, oI.raan, oI.argp, oI.nu, mu);
r1.initial = oI;
r1.final = parkingOrbit(cfg, vInf, oI);
[r1.cases, r1.table] = evaluateAll(cfg, oI, r1.final, mu);
r1.selected = selectCase(r1.table, cfg.sc1.selection);
r1.best = r1.cases(r1.selected);
r1.best.time = min(r1.best.dtTotal);
end

function park = parkingOrbit(cfg, vInf, oI)
p = cfg.sc1.parking;
switch p.mode
    case 'designed'
        park = designParkingOrbit(vInf, cfg.const.rEarth + p.periapsisAltitude, ...
            cfg.const.rEarth + p.apoapsisAltitude, oI, cfg.const.muEarth);
    case 'given'
        park = car2kep(p.state.r, p.state.v, cfg.const.muEarth);
    otherwise
        error('runScenario1:parkingMode', 'cfg.sc1.parking.mode must be ''designed'' or ''given''.');
end
end

function [cases, T] = evaluateAll(cfg, oI, oF, mu)
cases = [];
for s = 1:numel(cfg.sc1.strategies)
    for t = 1:numel(cfg.sc1.types)
        cases = [cases; evaluateStrategy(cfg.sc1.strategies{s}, cfg.sc1.types{t}, oI, oF, mu)]; %#ok<AGROW>
    end
end
dtHours = reshape([cases.dtTotal], 2, []).'/3600;
T = table(string({cases.name}).', [cases.dvTotal].', min(dtHours, [], 2), ...
    dtHours(:, 1), dtHours(:, 2), [cases.feasible].', ...
    'VariableNames', {'Strategy', 'DeltaV_kms', 'Time_h', 'TimeOption1_h', 'TimeOption2_h', 'Feasible'});
end

function k = selectCase(T, selection)
switch selection
    case 'minDeltaV'
        [~, k] = min(T.DeltaV_kms);
    case 'minTime'
        [~, k] = min(T.Time_h);
    otherwise
        k = find(T.Strategy == selection, 1);
        if isempty(k)
            error('runScenario1:badSelection', 'Unknown selection "%s".', selection);
        end
end
end
