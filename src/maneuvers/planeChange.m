function [dv, argpF, nu] = planeChange(a, e, iI, raanI, argpI, iF, raanF, mu)
%PLANECHANGE Single-impulse change of orbital plane (spherical triangle).
%   [dv, argpF, nu] = PLANECHANGE(a, e, iI, raanI, argpI, iF, raanF, mu)
%   rotates the orbit (a, e) from the plane (iI, raanI) to (iF, raanF) at
%   the intersection of the two planes. Of the two nodes, the one farther
%   from periapsis (lower speed, cheaper) is used. Returns the impulse dv
%   [km/s], the new argument of periapsis argpF and the true anomaly nu of
%   the maneuver [rad].

dRaan = mod(raanF - raanI + pi, 2*pi) - pi;     % in [-pi, pi)
di = iF - iI;

alpha = acos(clamp(cos(iI)*cos(iF) + sin(iI)*sin(iF)*cos(dRaan)));   % angle between planes
if alpha < 1e-12
    dv = 0; argpF = argpI; nu = 0;
    return
end
sa = sin(alpha);
if abs(sa) < 1e-12 || abs(sin(iI)) < 1e-12 || abs(sin(iF)) < 1e-12
    error('planeChange:degenerate', 'Degenerate geometry (alpha = 0/pi or equatorial orbit).');
end

% Argument of latitude of the node on the initial (uI) and final (uF) plane
sUI = clamp(sin(abs(dRaan))*sin(iF)/sa);
sUF = clamp(sin(abs(dRaan))*sin(iI)/sa);
if di >= 0
    cUI = clamp((-cos(iF) + cos(alpha)*cos(iI))/(sa*sin(iI)));
    cUF = clamp(( cos(iI) - cos(alpha)*cos(iF))/(sa*sin(iF)));
else
    cUI = clamp(( cos(iF) - cos(alpha)*cos(iI))/(sa*sin(iI)));
    cUF = clamp((-cos(iI) + cos(alpha)*cos(iF))/(sa*sin(iF)));
end
uI = mod(atan2(sUI, cUI), 2*pi);
uF = mod(atan2(sUF, cUF), 2*pi);

if (dRaan >= 0) == (di >= 0)
    nu    = mod(uI - argpI, 2*pi);
    argpF = mod(uF - nu, 2*pi);
else
    nu    = mod(2*pi - uI - argpI, 2*pi);
    argpF = mod(2*pi - uF - nu, 2*pi);
end

% Use the other node if this one is on the periapsis side
if cos(nu) > 0
    nu = mod(nu + pi, 2*pi);
end

p = a*(1 - e^2);
vTheta = sqrt(mu/p)*(1 + e*cos(nu));
dv = 2*vTheta*sin(alpha/2);
end

function x = clamp(x)
x = max(-1, min(1, x));
end
