function leg = hyperbolicLeg(vInf, rp, mu, rBody, rSoi)
%HYPERBOLICLEG Planar hyperbola with excess speed vInf and periapsis rp.
%   leg = HYPERBOLICLEG(vInf, rp, mu, rBody, rSoi) returns semi-major axis
%   a (< 0), eccentricity e, semi-latus rectum p, periapsis speed vp, local
%   circular and escape speeds, impact parameter, turn angle delta and the
%   true anomaly of the asymptote. leg.feasible is false if rp is below the
%   body surface rBody or outside its sphere of influence rSoi.

if nargin < 4 || isempty(rBody), rBody = 0;   end
if nargin < 5 || isempty(rSoi),  rSoi  = Inf; end
leg.vInf = vInf;
leg.rp = rp;
leg.a = -mu/vInf^2;
leg.e = 1 - rp/leg.a;
leg.p = leg.a*(1 - leg.e^2);
leg.vEsc = sqrt(2*mu/rp);
leg.vCirc = sqrt(mu/rp);
leg.vp = sqrt(vInf^2 + leg.vEsc^2);
leg.impactParameter = -leg.a*sqrt(leg.e^2 - 1);
leg.turnAngle = 2*asin(1/leg.e);
leg.nuInf = acos(-1/leg.e);
leg.feasible = true;
leg.reason = "";
if ~(isfinite(vInf) && vInf > 0)
    leg.feasible = false; leg.reason = "invalid v_inf";
elseif rp <= rBody
    leg.feasible = false; leg.reason = "periapsis below the surface";
elseif rp >= rSoi
    leg.feasible = false; leg.reason = "periapsis outside the SOI";
end
end
