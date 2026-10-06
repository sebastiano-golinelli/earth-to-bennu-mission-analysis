function [dv, nuI, nuF] = periapsisChange(a, e, argpI, argpF, mu)
%PERIAPSISCHANGE Single-impulse rotation of the apse line in the orbit plane.
%   [dv, nuI, nuF] = PERIAPSISCHANGE(a, e, argpI, argpF, mu) rotates the
%   argument of periapsis from argpI to argpF keeping a and e. The maneuver
%   can be done at two points: nuI and nuF (2x1) are their true anomalies on
%   the initial and on the final orbit [rad]. dv [km/s] is the same for both.

if e < 1e-10
    error('periapsisChange:circular', 'Argument of periapsis undefined for a circular orbit.');
end
if e >= 1
    error('periapsisChange:notElliptic', 'Only elliptic orbits (e < 1) are supported.');
end
dArgp = mod(argpF - argpI, 2*pi);
if dArgp < 1e-12 || abs(dArgp - 2*pi) < 1e-12
    dv = 0; nuI = [0; pi]; nuF = [0; pi];
    return
end
nuI = mod([dArgp/2; pi + dArgp/2], 2*pi);
nuF = mod([2*pi - dArgp/2; pi - dArgp/2], 2*pi);
p = a*(1 - e^2);
dv = 2*e*sqrt(mu/p)*abs(sin(dArgp/2));
end
