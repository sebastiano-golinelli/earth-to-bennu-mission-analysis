function c = evaluateStrategy(strategy, type, oI, oF, mu)
%EVALUATESTRATEGY Cost and duration of one geocentric transfer strategy.
%   c = EVALUATESTRATEGY(strategy, type, oI, oF, mu) transfers from the
%   initial orbit oI (struct a, e, i, raan, argp, nu = current position) to
%   the final orbit oF (same fields, nu = target position) using
%     strategy 'SPW'   : Shape (bitangent) -> Plane change -> periapsis (W)
%              'PSW'   : Plane change -> Shape -> periapsis
%              'PWS'   : Plane change -> periapsis -> Shape
%              'S(P)W' : Shape with the plane change inside the transfer arc -> periapsis
%     type     'ap', 'pa', 'pp', 'aa' : apses used by the bitangent transfer.
%
%   c.dv      = [dv1 bitangent, dv2 bitangent, dv plane, dv periapsis] [km/s]
%   c.dvTotal = sum of |dv| [km/s]
%   c.dtTotal = [1x2] total time [s] for the two possible periapsis-change
%               points (coasting phases included)
%   c.geo     = intermediate arguments of periapsis, used for plotting

c = struct('name', [strategy ' ' type], 'strategy', strategy, 'type', type, ...
    'feasible', true, 'dv', nan(1, 4), 'dvTotal', NaN, 'dtTotal', nan(1, 2), 'geo', struct());

nuDep = apsisAnomaly(type(1));     % departure apsis on the initial orbit
nuArr = apsisAnomaly(type(2));     % arrival apsis on the final orbit
tof = @(a, e, nu1, nu2) timeOfFlight(a, e, nu1, nu2, mu);

switch strategy
    case 'SPW'
        dt1 = tof(oI.a, oI.e, oI.nu, nuDep);
        [dv1, dv2, dt2, argpS] = bitangentTransfer(oI.a, oI.e, oF.a, oF.e, type, oI.argp, mu);
        [dvP, argpP, nuP] = planeChange(oF.a, oF.e, oI.i, oI.raan, argpS, oF.i, oF.raan, mu);
        dt3 = tof(oF.a, oF.e, nuArr, nuP);
        [dvW, nuWi, nuWf] = periapsisChange(oF.a, oF.e, argpP, oF.argp, mu);
        dt = dt1 + dt2 + dt3 + tof(oF.a, oF.e, nuP, nuWi) + tof(oF.a, oF.e, nuWf, oF.nu);
        c.geo = struct('argpShape', argpS, 'argpPlane', argpP);

    case 'PSW'
        [dvP, argpP, nuP] = planeChange(oI.a, oI.e, oI.i, oI.raan, oI.argp, oF.i, oF.raan, mu);
        dt1 = tof(oI.a, oI.e, oI.nu, nuP);
        dt2 = tof(oI.a, oI.e, nuP, nuDep);
        [dv1, dv2, dt3, argpS] = bitangentTransfer(oI.a, oI.e, oF.a, oF.e, type, argpP, mu);
        [dvW, nuWi, nuWf] = periapsisChange(oF.a, oF.e, argpS, oF.argp, mu);
        dt = dt1 + dt2 + dt3 + tof(oF.a, oF.e, nuArr, nuWi) + tof(oF.a, oF.e, nuWf, oF.nu);
        c.geo = struct('argpPlane', argpP, 'argpShape', argpS);

    case 'PWS'
        [dvP, argpP, nuP] = planeChange(oI.a, oI.e, oI.i, oI.raan, oI.argp, oF.i, oF.raan, mu);
        % periapsis target before the bitangent (pp and aa rotate the apse line by pi)
        argpW = oF.argp;
        if type(1) == type(2), argpW = mod(oF.argp - pi, 2*pi); end
        dt1 = tof(oI.a, oI.e, oI.nu, nuP);
        [dvW, nuWi, nuWf] = periapsisChange(oI.a, oI.e, argpP, argpW, mu);
        dt2 = tof(oI.a, oI.e, nuP, nuWi);
        [dv1, dv2, dt3] = bitangentTransfer(oI.a, oI.e, oF.a, oF.e, type, argpW, mu);
        dt4 = tof(oI.a, oI.e, nuWf, nuDep);
        dt5 = tof(oF.a, oF.e, nuArr, oF.nu);
        dt = dt1 + dt2 + dt4 + dt3 + dt5;
        c.geo = struct('argpPlane', argpP, 'argpPeriapsis', argpW);

    case 'S(P)W'
        dt1 = tof(oI.a, oI.e, oI.nu, nuDep);
        m = bitangentWithPlaneChange(oI, oF, type, mu);
        if ~m.feasible
            c.feasible = false;
            return
        end
        dv1 = m.dv1; dv2 = m.dv2; dvP = m.dvPlane;
        [dvW, nuWi, nuWf] = periapsisChange(oF.a, oF.e, m.argpFinal, oF.argp, mu);
        dt = dt1 + sum(m.dtCoast) + tof(oF.a, oF.e, nuArr, nuWi) + tof(oF.a, oF.e, nuWf, oF.nu);
        c.geo = m;

    otherwise
        error('evaluateStrategy:badStrategy', 'Unknown strategy %s.', strategy);
end

c.dv = [dv1, dv2, dvP, dvW];
c.dvTotal = sum(abs(c.dv));
c.dtTotal = dt(:).';
end
