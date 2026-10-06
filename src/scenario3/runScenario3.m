function r3 = runScenario3(cfg, r1, r2)
%RUNSCENARIO3 Patched conics: Earth escape and capture at the target.
%   r3 = RUNSCENARIO3(cfg, r1, r2) uses the parking orbit reached in
%   Scenario 1 (r1.final) and the excess velocities of the best heliocentric
%   transfer of Scenario 2 (r2.best) to compute
%     - the spheres of influence of Earth and target,
%     - the coplanar escape (tangential burn at parking periapsis),
%     - the capture into circular orbits at several altitudes,
%     - the non-coplanar escape, with the burn point optimized along the
%       parking orbit (excess velocity rotated from ecliptic to equator).

c = cfg.const;
s = cfg.sc3;
tgt = cfg.sc2.target;
park = r1.final;
park.rp = park.a*(1 - park.e);
park.ra = park.a*(1 + park.e);
park.vp = sqrt(c.muEarth*(2/park.rp - 1/park.a));
r3.parking = park;

best = r2.best;
r3.vInfDepEcliptic = best.vT1 - best.vEarth;
r3.vInfArr = best.vT2 - best.vTarget;
r3.vInfDep = norm(r3.vInfDepEcliptic);
r3.vInfArrMag = norm(r3.vInfArr);

% ---- Spheres of influence ---------------------------------------------------
r3.soiEarth  = sphereOfInfluence(cfg.sc2.earth.a, c.muEarth, c.muSun);
r3.soiTarget = sphereOfInfluence(tgt.a, tgt.mu, c.muSun);
r3.soiEarthAtPatch  = sphereOfInfluence(norm(best.r1), c.muEarth, c.muSun);
r3.soiTargetAtPatch = sphereOfInfluence(norm(best.r2), tgt.mu, c.muSun);
if strcmp(s.soiMode, 'patch')
    soiE = r3.soiEarthAtPatch; soiT = r3.soiTargetAtPatch;
else
    soiE = r3.soiEarth; soiT = r3.soiTarget;
end

% ---- Coplanar escape: tangential burn at parking periapsis ------------------
dep = hyperbolicLeg(r3.vInfDep, park.rp, c.muEarth, c.rEarth, soiE);
if ~dep.feasible
    error('runScenario3:escape', 'Escape hyperbola not feasible: %s.', dep.reason);
end
dep.dv = dep.vp - park.vp;
dep.timeInSoi = hyperbolicTimeOfFlight(dep.a, dep.e, c.muEarth, soiE);
r3.escape = dep;

% Sensitivity of the escape cost to the periapsis altitude (same parking a)
rpSweep = c.rEarth + s.escapeAltitudeSweep(:);
dvSweep = nan(size(rpSweep));
for k = 1:numel(rpSweep)
    leg = hyperbolicLeg(r3.vInfDep, rpSweep(k), c.muEarth, c.rEarth, soiE);
    dvSweep(k) = leg.vp - sqrt(c.muEarth*(2/rpSweep(k) - 1/park.a));
end
r3.escapeSweep = table(s.escapeAltitudeSweep(:), rpSweep, dvSweep, ...
    'VariableNames', {'Altitude_km', 'Periapsis_km', 'DeltaV_kms'});

% ---- Capture at the target --------------------------------------------------
h = s.captureAltitudes(:);
nCap = numel(h);
[dvCirc, dvMarg, ecc, ratio, tide, srp] = deal(nan(nCap, 1));
level = strings(nCap, 1);
for k = 1:nCap
    rCap = tgt.radius + h(k);
    leg = hyperbolicLeg(r3.vInfArrMag, rCap, tgt.mu, tgt.radius, soiT);
    [~, info] = captureCost(r3.vInfArrMag, rCap, tgt.mu);
    v = captureValidity(rCap, soiT, tgt.mu, cfg, norm(best.r2));
    dvCirc(k) = info.dvCircular; dvMarg(k) = info.dvMarginal; ecc(k) = leg.e;
    ratio(k) = v.soiRatio; tide(k) = v.tideRatio; srp(k) = v.srpRatio; level(k) = v.level;
end
r3.captureTable = table(h, tgt.radius + h, ecc, dvCirc, dvMarg, ratio, level, tide, srp, ...
    'VariableNames', {'Altitude_km', 'Radius_km', 'HyperbolaEcc', 'DeltaVCircular_kms', ...
    'DeltaVMarginal_kms', 'RadiusOverSOI', 'PatchValidity', 'TideRatio', 'SRPRatio'});
[~, kNom] = min(abs(h - s.nominalCaptureAltitude));
r3.capture = table2struct(r3.captureTable(kNom, :));
legNom = hyperbolicLeg(r3.vInfArrMag, tgt.radius + h(kNom), tgt.mu);
r3.capture.timeInSoi = hyperbolicTimeOfFlight(legNom.a, legNom.e, tgt.mu, soiT);
r3.capture.dv = r3.capture.DeltaVCircular_kms;

% ---- Non-coplanar escape ----------------------------------------------------
r3.vInfDepEquatorial = eclipticToEquatorial(r3.vInfDepEcliptic, c.obliquity);

nuBurn = linspace(0, 2*pi, s.burnSweepPoints);
dvBurn = nan(size(nuBurn));
for k = 1:numel(nuBurn)
    [rb, vb] = kep2car(park.a, park.e, park.i, park.raan, park.argp, nuBurn(k), c.muEarth);
    H = escapeHyperbola3D(r3.vInfDepEquatorial, rb, vb, c.muEarth, c.rEarth);
    if H.feasible, dvBurn(k) = H.dv; end
end
if all(isnan(dvBurn))
    error('runScenario3:noncoplanar', 'No feasible non-coplanar escape along the parking orbit.');
end
[~, kOpt] = min(dvBurn);
[rb, vb] = kep2car(park.a, park.e, park.i, park.raan, park.argp, nuBurn(kOpt), c.muEarth);
r3.escape3D = escapeHyperbola3D(r3.vInfDepEquatorial, rb, vb, c.muEarth, c.rEarth);
r3.escape3D.nuBurnOnParking = nuBurn(kOpt);
r3.burnSweep = struct('nu', nuBurn, 'dv', dvBurn);

% ---- Totals -----------------------------------------------------------------
r3.dvCoplanar = r3.escape.dv + r3.capture.dv;
r3.dvNonCoplanar = r3.escape3D.dv + r3.capture.dv;
end
