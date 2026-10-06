function [dv, info] = captureCost(vInf, rp, mu, ra)
%CAPTURECOST Impulsive capture at the periapsis of the arrival hyperbola.
%   [dv, info] = CAPTURECOST(vInf, rp, mu, ra) inserts into an orbit with
%   periapsis rp and apoapsis ra (default ra = rp: circular orbit).
%   info.dvCircular and info.dvMarginal are the costs of a circular capture
%   and of a barely closed orbit (ra -> infinity).

if nargin < 4 || isempty(ra), ra = rp; end
vpHyp = sqrt(vInf^2 + 2*mu/rp);
vpOrbit = sqrt(mu*(2/rp - 2/(rp + ra)));
dv = abs(vpHyp - vpOrbit);
info.vpHyperbola = vpHyp;
info.dvCircular = vpHyp - sqrt(mu/rp);
info.dvMarginal = vpHyp - sqrt(2*mu/rp);
end
