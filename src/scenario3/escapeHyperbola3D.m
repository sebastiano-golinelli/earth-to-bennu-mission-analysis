function H = escapeHyperbola3D(vInf, rBurn, vPark, mu, rBody)
%ESCAPEHYPERBOLA3D Non-coplanar escape hyperbola from a given burn point.
%   H = ESCAPEHYPERBOLA3D(vInf, rBurn, vPark, mu, rBody) finds the escape
%   hyperbola that passes through the burn point rBurn (on the parking
%   orbit, where the velocity is vPark) and whose outgoing asymptote is
%   parallel to the required excess velocity vInf (all 3x1, same inertial
%   frame). The periapsis of the hyperbola is generally not at the burn
%   point. H.dv = |v_hyperbola(rBurn) - vPark|. H.feasible is false if no
%   such hyperbola exists or its periapsis is below the surface rBody.

if nargin < 5 || isempty(rBody), rBody = 0; end
vInf = vInf(:); rBurn = rBurn(:); vPark = vPark(:);
H = struct('a', NaN, 'e', NaN, 'i', NaN, 'raan', NaN, 'argp', NaN, 'nuBurn', NaN, ...
    'nuInf', NaN, 'rp', NaN, 'alpha', NaN, 'rBurn', rBurn, 'vBurn', nan(3, 1), ...
    'dvVec', nan(3, 1), 'dv', NaN, 'feasible', false, 'reason', "");

vInfMag = norm(vInf);
rB = norm(rBurn);
if vInfMag <= 0 || rB <= rBody
    H.reason = "invalid input"; return
end
uInf = vInf/vInfMag;
uR = rBurn/rB;
a = -mu/vInfMag^2;
alpha = acos(max(-1, min(1, dot(uR, uInf))));    % angle between burn point and asymptote
H.a = a; H.alpha = alpha;

% Eccentricity from rB = a(1-e^2)/(1 + e cos(nuBurn)), nuBurn = acos(-1/e) - alpha
eMax = 1 - rB/a;                                  % periapsis cannot exceed rB
if eMax <= 1.000001
    H.reason = "no admissible eccentricity"; return
end
f = @(e) rB - a.*(1 - e.^2)./(1 + e.*cos(acos(-1./e) - alpha));
eGrid = linspace(1 + 1e-6, eMax - 1e-9, 300);
fGrid = f(eGrid);
k = find(fGrid(1:end-1).*fGrid(2:end) < 0, 1);
if isempty(k)
    H.reason = "incompatible geometry"; return
end
e = fzero(f, eGrid([k k+1]));
H.e = e;
H.nuInf = acos(-1/e);
H.nuBurn = H.nuInf - alpha;
H.rp = a*(1 - e);

% Hyperbola plane = span(rBurn, vInf)
h = cross(uR, uInf);
if norm(h) < 1e-8
    H.reason = "burn point aligned with v_inf"; return
end
h = h/norm(h);
H.i = acos(max(-1, min(1, h(3))));
n = cross([0; 0; 1], h);
if norm(n) < 1e-12
    n = [1; 0; 0]; H.raan = 0;
else
    n = n/norm(n); H.raan = mod(atan2(n(2), n(1)), 2*pi);
end
u = atan2(dot(cross(n, uR), h), dot(n, uR));     % argument of latitude of the burn point
H.argp = mod(u - H.nuBurn, 2*pi);

[rCheck, H.vBurn] = kep2car(a, e, H.i, H.raan, H.argp, H.nuBurn, mu);
if norm(rCheck - rBurn) > 1e-6*rB
    H.reason = "geometric inconsistency"; return
end
H.dvVec = H.vBurn - vPark;
H.dv = norm(H.dvVec);
H.feasible = H.rp > rBody;
if ~H.feasible
    H.reason = "periapsis below the surface";
end
end
