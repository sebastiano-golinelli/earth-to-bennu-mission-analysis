function t = hyperbolicTimeOfFlight(a, e, mu, r)
%HYPERBOLICTIMEOFFLIGHT Time [s] from periapsis to radius r on a hyperbola (a < 0).

F = acosh(max(1, (1 - r/a)/e));          % r = a(1 - e cosh F)
t = (e*sinh(F) - F)/sqrt(mu/(-a)^3);
end
