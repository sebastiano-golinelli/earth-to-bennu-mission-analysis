function P = conicPoints(a, e, i, raan, argp, nu)
%CONICPOINTS Inertial positions along a conic for an array of true anomalies.
%   P = CONICPOINTS(a, e, i, raan, argp, nu) returns a 3xN matrix of
%   positions [km]. Used for plotting; works for ellipses and hyperbolas.

nu = nu(:).';
p = a*(1 - e^2);
r = p./(1 + e*cos(nu));
P = perifocalToInertial(raan, i, argp) * [r.*cos(nu); r.*sin(nu); zeros(size(nu))];
end
