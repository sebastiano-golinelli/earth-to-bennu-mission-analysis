function plotStrategySequence(c, oI, oF, dnu)
%PLOTSTRATEGYSEQUENCE Orbits visited by one geocentric transfer strategy.
%   PLOTSTRATEGYSEQUENCE(c, oI, oF, dnu) draws on the current axes the
%   launch orbit oI, the intermediate orbits and transfer arcs of strategy
%   c (output of evaluateStrategy) and the parking orbit oF.

g = c.geo;
col = lines(5);
plotOrbitArc(oI.a, oI.e, oI.i, oI.raan, oI.argp, 0, 2*pi, dnu, col(1, :), 'Launch orbit (GTO)', '-', 1.5);
switch c.strategy
    case 'SPW'
        bitangentArc(oI.a, oI.e, oI.i, oI.raan, oI.argp, oF.a, oF.e, c.type, dnu, col(2, :));
        plotOrbitArc(oF.a, oF.e, oI.i, oI.raan, g.argpShape, 0, 2*pi, dnu, col(3, :), 'After shape change', '-.', 1.3);
        plotOrbitArc(oF.a, oF.e, oF.i, oF.raan, g.argpPlane, 0, 2*pi, dnu, col(4, :), 'After plane change', ':', 1.6);
    case 'PSW'
        plotOrbitArc(oI.a, oI.e, oF.i, oF.raan, g.argpPlane, 0, 2*pi, dnu, col(2, :), 'After plane change', '-.', 1.3);
        bitangentArc(oI.a, oI.e, oF.i, oF.raan, g.argpPlane, oF.a, oF.e, c.type, dnu, col(3, :));
        plotOrbitArc(oF.a, oF.e, oF.i, oF.raan, g.argpShape, 0, 2*pi, dnu, col(4, :), 'After shape change', ':', 1.6);
    case 'PWS'
        plotOrbitArc(oI.a, oI.e, oF.i, oF.raan, g.argpPlane, 0, 2*pi, dnu, col(2, :), 'After plane change', '-.', 1.3);
        plotOrbitArc(oI.a, oI.e, oF.i, oF.raan, g.argpPeriapsis, 0, 2*pi, dnu, col(3, :), 'After periapsis change', ':', 1.6);
        bitangentArc(oI.a, oI.e, oF.i, oF.raan, g.argpPeriapsis, oF.a, oF.e, c.type, dnu, col(4, :));
    case 'S(P)W'
        if strcmpi(c.type, 'ap'), nuDep = pi; nuArr = 2*pi; else, nuDep = 0; nuArr = pi; end
        plotOrbitArc(g.aT, g.eT, oI.i, oI.raan, g.argpBefore, nuDep, g.nuNode, dnu, col(2, :), ...
            'Transfer arc (launch plane)', '--', 1.8);
        plotOrbitArc(g.aT, g.eT, oF.i, oF.raan, g.argpAfter, g.nuNode, nuArr, dnu, col(3, :), ...
            'Transfer arc (parking plane)', '--', 1.8);
        plotOrbitArc(oF.a, oF.e, oF.i, oF.raan, g.argpFinal, 0, 2*pi, dnu, col(4, :), 'After transfer', ':', 1.6);
        P = conicPoints(g.aT, g.eT, oI.i, oI.raan, g.argpBefore, g.nuNode);
        plot3(P(1), P(2), P(3), 'p', 'Color', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', 11, ...
            'DisplayName', 'Plane change');
end
plotOrbitArc(oF.a, oF.e, oF.i, oF.raan, oF.argp, 0, 2*pi, dnu, col(5, :), 'Parking orbit', '-', 2.2);
end

function bitangentArc(a1, e1, i, raan, argp1, a2, e2, type, dnu, color)
% Half ellipse between the departure apsis of orbit 1 and the arrival apsis
% of orbit 2 (same plane, apse line rotated by pi for pp and aa).
nuDep = apsisAnomaly(type(1));
nuArr = apsisAnomaly(type(2));
rDep = a1*(1 - e1^2)/(1 + e1*cos(nuDep));
rArr = a2*(1 - e2^2)/(1 + e2*cos(nuArr));
aT = (rDep + rArr)/2;
eT = abs(rArr - rDep)/(rArr + rDep);
if rDep <= rArr
    plotOrbitArc(aT, eT, i, raan, mod(argp1 + nuDep, 2*pi), 0, pi, dnu, color, 'Bitangent transfer', '--', 1.8);
else
    plotOrbitArc(aT, eT, i, raan, mod(argp1 + nuDep - pi, 2*pi), pi, 2*pi, dnu, color, 'Bitangent transfer', '--', 1.8);
end
end
