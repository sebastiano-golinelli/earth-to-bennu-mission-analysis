function g = transferGeometry(x, cfg)
%TRANSFERGEOMETRY Conic through the departure and arrival points.
%   g = TRANSFERGEOMETRY(x, cfg) with x = [nuDep, nuArr, argpT] [rad]:
%   departure true anomaly on the Earth orbit, arrival true anomaly on the
%   target orbit, argument of periapsis of the transfer orbit.
%   The transfer plane contains r1 and r2. With cfg.sc2.allowLongArc = true
%   the motion is always prograde, so the arc from r1 to r2 can be shorter
%   or longer than 180 deg; with false the angular momentum is r1 x r2
%   (short arc only, possibly retrograde).
%   g.ok = false when the geometry is degenerate (r1 and r2 aligned).

mu = cfg.const.muSun;
oE = cfg.sc2.earth;
oT = cfg.sc2.target;
g.x = mod(x(:).', 2*pi);
g.ok = false;
[g.r1, g.vEarth]  = kep2car(oE.a, oE.e, oE.i, oE.raan, oE.argp, g.x(1), mu);
[g.r2, g.vTarget] = kep2car(oT.a, oT.e, oT.i, oT.raan, oT.argp, g.x(2), mu);
r1n = norm(g.r1);
r2n = norm(g.r2);

h = cross(g.r1, g.r2);
if norm(h) < 1e-8*r1n*r2n
    return
end
h = h/norm(h);
if cfg.sc2.allowLongArc && h(3) < 0
    h = -h;
end
g.iT = acos(max(-1, min(1, h(3))));
n = cross([0; 0; 1], h);
if norm(n) < 1e-12
    g.raanT = 0;
else
    n = n/norm(n);
    g.raanT = mod(atan2(n(2), n(1)), 2*pi);
end

% True anomalies of r1 and r2 on the transfer orbit
Q = perifocalToInertial(g.raanT, g.iT, g.x(3));
r1pf = Q.'*g.r1;
r2pf = Q.'*g.r2;
g.nu1T = mod(atan2(r1pf(2), r1pf(1)), 2*pi);
g.nu2T = mod(atan2(r2pf(2), r2pf(1)), 2*pi);

% Conic r = p/(1 + e cos(nu)) written at r1 and r2
c1 = cos(g.nu1T);
c2 = cos(g.nu2T);
den = r1n*c1 - r2n*c2;
if abs(den) < 1e-10*max(r1n, r2n)
    return
end
g.eT = (r2n - r1n)/den;
g.pT = r1n*(1 + g.eT*c1);
g.ok = true;
end
