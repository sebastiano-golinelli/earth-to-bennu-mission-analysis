function park = designParkingOrbit(vInf, rp, ra, ref, mu)
%DESIGNPARKINGORBIT Parking orbit aligned with the departure asymptote.
%   park = DESIGNPARKINGORBIT(vInf, rp, ra, ref, mu) returns the elements
%   (struct a, e, i, raan, argp, nu) of an orbit with periapsis rp and
%   apoapsis ra [km] such that a tangential burn at periapsis puts the
%   spacecraft on an escape hyperbola whose outgoing asymptote is parallel
%   to vInf (3x1, equatorial frame, km/s):
%     - the orbit plane contains vInf; the inclination is the one of the
%       reference orbit ref (e.g. the launch orbit) if possible, otherwise
%       the minimum one, equal to the declination of vInf;
%     - of the two planes with that inclination, the one closer to the
%       plane of ref is chosen (cheaper plane change);
%     - the periapsis precedes the asymptote by the true anomaly of the
%       asymptote of the escape hyperbola.
%   nu = 0: the spacecraft must arrive at periapsis, where the escape burn
%   is performed.

u = vInf(:)/norm(vInf);
dec = asin(u(3));                  % declination of the asymptote
ra0 = atan2(u(2), u(1));           % right ascension of the asymptote
i = max(ref.i, abs(dec));
k = max(-1, min(1, tan(dec)/tan(i)));
raanCand = mod([ra0 - asin(k), ra0 - pi + asin(k)], 2*pi);
planeAngle = acos(max(-1, min(1, cos(ref.i)*cos(i) + sin(ref.i)*sin(i)*cos(raanCand - ref.raan))));
[~, kBest] = min(planeAngle);
raan = raanCand(kBest);

% Argument of latitude of the asymptote in the parking plane
n = [cos(raan); sin(raan); 0];
h = [sin(raan)*sin(i); -cos(raan)*sin(i); cos(i)];
uInf = atan2(dot(u, cross(h, n)), dot(u, n));

eHyp = 1 + rp*norm(vInf)^2/mu;
nuInf = acos(-1/eHyp);
park = struct('a', (rp + ra)/2, 'e', (ra - rp)/(ra + rp), 'i', i, 'raan', raan, ...
    'argp', mod(uInf - nuInf, 2*pi), 'nu', 0);
park.asymptote = struct('rightAscension', ra0, 'declination', dec);
end
