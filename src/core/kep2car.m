function [r, v] = kep2car(a, e, i, raan, argp, nu, mu)
%KEP2CAR Keplerian elements to inertial position and velocity.
%   [r, v] = KEP2CAR(a, e, i, raan, argp, nu, mu) returns the position r [km]
%   and velocity v [km/s] (3x1) on a conic with semi-major axis a [km]
%   (negative for hyperbolas), eccentricity e, inclination i, right
%   ascension of the ascending node raan, argument of periapsis argp and
%   true anomaly nu [rad]. mu is the gravitational parameter [km^3/s^2].

p = a*(1 - e^2);
rPf = p/(1 + e*cos(nu)) * [cos(nu); sin(nu); 0];
vPf = sqrt(mu/p) * [-sin(nu); e + cos(nu); 0];
Q = perifocalToInertial(raan, i, argp);
r = Q*rPf;
v = Q*vPf;
end
