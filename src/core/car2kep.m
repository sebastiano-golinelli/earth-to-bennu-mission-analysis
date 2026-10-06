function kep = car2kep(r, v, mu)
%CAR2KEP Inertial position and velocity to Keplerian elements.
%   kep = CAR2KEP(r, v, mu) returns a struct with fields a [km], e, i, raan,
%   argp, nu [rad] for the state r [km], v [km/s] and gravitational
%   parameter mu [km^3/s^2]. Equatorial or circular orbits, where raan or
%   argp are undefined, raise an error.

r = r(:); v = v(:);
rn = norm(r);
h = cross(r, v);
n = cross([0; 0; 1], h);
eVec = cross(v, h)/mu - r/rn;
e = norm(eVec);
if norm(n) < 1e-12*norm(h) || e < 1e-12
    error('car2kep:singular', 'Equatorial or circular orbit: RAAN or argument of periapsis undefined.');
end
n = n/norm(n);

kep.a = 1/(2/rn - dot(v, v)/mu);
kep.e = e;
kep.i = acos(clamp(h(3)/norm(h)));
kep.raan = acos(clamp(n(1)));
if n(2) < 0, kep.raan = 2*pi - kep.raan; end
kep.argp = acos(clamp(dot(n, eVec)/e));
if eVec(3) < 0, kep.argp = 2*pi - kep.argp; end
kep.nu = acos(clamp(dot(r, eVec)/(rn*e)));
if dot(r, v) < 0, kep.nu = 2*pi - kep.nu; end
end

function x = clamp(x)
x = max(-1, min(1, x));
end
