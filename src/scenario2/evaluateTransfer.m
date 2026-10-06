function sol = evaluateTransfer(x, cfg)
%EVALUATETRANSFER Delta-v and time of flight of a two-impulse transfer.
%   sol = EVALUATETRANSFER(x, cfg) with x = [nuDep, nuArr, argpT] [rad]
%   (see transferGeometry). sol.valid is false, with a positive penalty,
%   when the transfer is not an ellipse or its perihelion is below
%   cfg.sc2.minPerihelion. The cost is J = dv + cfg.sc2.weightTOF*TOF[days].

sol = struct('x', [NaN NaN NaN], 'J', Inf, 'dv', NaN, 'dv1', NaN, 'dv2', NaN, ...
    'tofDays', NaN, 'aT', NaN, 'eT', NaN, 'iT', NaN, 'raanT', NaN, 'argpT', NaN, ...
    'nu1T', NaN, 'nu2T', NaN, 'arcDeg', NaN, 'rpT', NaN, 'valid', false, 'penalty', 0, ...
    'r1', nan(3, 1), 'r2', nan(3, 1), 'vEarth', nan(3, 1), 'vTarget', nan(3, 1), ...
    'vT1', nan(3, 1), 'vT2', nan(3, 1), 'method', "");

mu = cfg.const.muSun;
g = transferGeometry(x, cfg);
sol.x = g.x;
sol.argpT = g.x(3);
sol.r1 = g.r1; sol.r2 = g.r2; sol.vEarth = g.vEarth; sol.vTarget = g.vTarget;
if ~g.ok
    sol.penalty = 1e8;
    return
end
sol.iT = g.iT; sol.raanT = g.raanT; sol.nu1T = g.nu1T; sol.nu2T = g.nu2T;
sol.arcDeg = rad2deg(mod(g.nu2T - g.nu1T, 2*pi));
eT = g.eT; pT = g.pT;
if ~isfinite(eT) || ~isfinite(pT)
    sol.penalty = 1e9;
    return
end

penalty = 0;
if eT < 0,  penalty = penalty + 1e6 + 1e6*abs(eT);   end
if eT >= 1, penalty = penalty + 1e6 + 1e6*(eT - 1)^2; end
if pT <= 0, penalty = penalty + 1e8 + 1e-2*abs(pT);  end
if penalty > 0
    sol.eT = eT;
    sol.penalty = penalty;
    return
end

aT = pT/(1 - eT^2);
sol.aT = aT; sol.eT = eT; sol.rpT = aT*(1 - eT);
if sol.rpT <= cfg.sc2.minPerihelion
    sol.penalty = 1e7 + (cfg.sc2.minPerihelion - sol.rpT)^2;
    return
end

[~, sol.vT1] = kep2car(aT, eT, g.iT, g.raanT, g.x(3), g.nu1T, mu);
[~, sol.vT2] = kep2car(aT, eT, g.iT, g.raanT, g.x(3), g.nu2T, mu);
sol.dv1 = norm(sol.vT1 - g.vEarth);
sol.dv2 = norm(g.vTarget - sol.vT2);
sol.dv = sol.dv1 + sol.dv2;
sol.tofDays = timeOfFlight(aT, eT, g.nu1T, g.nu2T, mu)/86400;
sol.J = sol.dv + cfg.sc2.weightTOF*sol.tofDays;
sol.valid = isfinite(sol.J) && sol.tofDays > 0;
end
