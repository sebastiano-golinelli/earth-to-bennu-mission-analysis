function dt = timeOfFlight(a, e, nu1, nu2, mu)
%TIMEOFFLIGHT Time of flight between two true anomalies on an ellipse.
%   dt = TIMEOFFLIGHT(a, e, nu1, nu2, mu) returns the time [s] needed to go
%   from true anomaly nu1 to nu2 [rad] in the direction of motion, so that
%   0 <= dt < orbital period (dt = 0 when nu1 == nu2). nu1 and nu2 can be
%   arrays of compatible size.

if e >= 1
    error('timeOfFlight:notElliptic', 'Only elliptic orbits (e < 1) are supported.');
end
nu1 = mod(nu1, 2*pi);
nu2 = mod(nu2, 2*pi);
E1 = mod(2*atan2(sqrt(1 - e)*sin(nu1/2), sqrt(1 + e)*cos(nu1/2)), 2*pi);
E2 = mod(2*atan2(sqrt(1 - e)*sin(nu2/2), sqrt(1 + e)*cos(nu2/2)), 2*pi);
M1 = E1 - e*sin(E1);
M2 = E2 - e*sin(E2);
dt = mod(M2 - M1, 2*pi)/sqrt(mu/a^3);
end
