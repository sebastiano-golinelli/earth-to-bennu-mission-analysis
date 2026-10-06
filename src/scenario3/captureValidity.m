function v = captureValidity(rp, rSoi, muBody, cfg, distanceToSun)
%CAPTUREVALIDITY How well the two-body capture model holds at radius rp.
%   v = CAPTUREVALIDITY(rp, rSoi, muBody, cfg, distanceToSun) returns
%     v.soiRatio   rp/rSoi (patched-conic validity: lower is better)
%     v.level      "high" (< 1/3), "medium" (< 1/2), "low" (< 1), "outside SOI"
%     v.tideRatio  solar tidal acceleration / body gravity at rp
%     v.srpRatio   solar radiation pressure acceleration / body gravity at rp,
%                  for the spacecraft defined in cfg.sc3.srp
%   Ratios above ~0.1 mean the perturbation is not negligible.

v.soiRatio = rp/rSoi;
if v.soiRatio >= 1
    v.level = "outside SOI";
elseif v.soiRatio > 0.5
    v.level = "low";
elseif v.soiRatio > 1/3
    v.level = "medium";
else
    v.level = "high";
end
gBody = muBody/rp^2;                                      % [km/s^2]
v.tideRatio = 2*cfg.const.muSun*rp/distanceToSun^3/gBody;
s = cfg.sc3.srp;
aSrp = s.Cr*s.pressure1AU*(cfg.const.AU/distanceToSun)^2*s.areaToMass/1000;   % [km/s^2]
v.srpRatio = aSrp/gBody;
end
