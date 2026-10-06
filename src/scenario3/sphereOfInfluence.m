function rSoi = sphereOfInfluence(distanceToSun, muBody, muSun)
%SPHEREOFINFLUENCE Laplace sphere of influence radius r = d*(muBody/muSun)^(2/5).
%   distanceToSun is the semi-major axis (classical definition) or the
%   actual distance at the patch point.

rSoi = distanceToSun*(muBody/muSun)^(2/5);
end
