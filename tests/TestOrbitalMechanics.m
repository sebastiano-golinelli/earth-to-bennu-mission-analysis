classdef TestOrbitalMechanics < matlab.unittest.TestCase
    %TESTORBITALMECHANICS Element conversions and time of flight.

    properties (Constant)
        Mu = 398600.4418
    end

    methods (TestClassSetup)
        function addPaths(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src'), 'IncludingSubfolders', true));
        end
    end

    methods (Test)
        function roundTripKeplerCartesian(tc)
            rng(1);
            for k = 1:200
                o = randomOrbit();
                [r, v] = kep2car(o.a, o.e, o.i, o.raan, o.argp, o.nu, tc.Mu);
                b = car2kep(r, v, tc.Mu);
                tc.verifyEqual(b.a, o.a, 'RelTol', 1e-10);
                tc.verifyEqual(b.e, o.e, 'AbsTol', 1e-10);
                tc.verifyEqual(b.i, o.i, 'AbsTol', 1e-9);
                tc.verifyEqual(angleDiff(b.raan, o.raan), 0, 'AbsTol', 1e-9);
                tc.verifyEqual(angleDiff(b.argp, o.argp), 0, 'AbsTol', 1e-8);
                tc.verifyEqual(angleDiff(b.nu, o.nu), 0, 'AbsTol', 1e-8);
            end
        end

        function energyAndAngularMomentum(tc)
            rng(2);
            for k = 1:50
                o = randomOrbit();
                [r, v] = kep2car(o.a, o.e, o.i, o.raan, o.argp, o.nu, tc.Mu);
                tc.verifyEqual(norm(v)^2/2 - tc.Mu/norm(r), -tc.Mu/(2*o.a), 'RelTol', 1e-12);
                tc.verifyEqual(norm(cross(r, v)), sqrt(tc.Mu*o.a*(1 - o.e^2)), 'RelTol', 1e-12);
            end
        end

        function timeOfFlightBasicProperties(tc)
            a = 24000; e = 0.7; T = 2*pi*sqrt(a^3/tc.Mu);
            tc.verifyEqual(timeOfFlight(a, e, 0, pi, tc.Mu), T/2, 'RelTol', 1e-12);
            tc.verifyEqual(timeOfFlight(a, e, 1.2, 1.2, tc.Mu), 0);
            t12 = timeOfFlight(a, e, 0.3, 2.5, tc.Mu);
            t23 = timeOfFlight(a, e, 2.5, 5.9, tc.Mu);
            tc.verifyEqual(t12 + t23, timeOfFlight(a, e, 0.3, 5.9, tc.Mu), 'RelTol', 1e-12);
            tc.verifyEqual(timeOfFlight(a, e, 5.9, 0.3, tc.Mu), T - timeOfFlight(a, e, 0.3, 5.9, tc.Mu), 'RelTol', 1e-12);
        end

        function timeOfFlightAgainstPropagation(tc)
            % Kepler's equation versus numerical integration of the two-body problem
            o = struct('a', 30000, 'e', 0.6, 'i', 0.5, 'raan', 1, 'argp', 2, 'nu', 0.4);
            nu2 = 4.0;
            dt = timeOfFlight(o.a, o.e, o.nu, nu2, tc.Mu);
            [r0, v0] = kep2car(o.a, o.e, o.i, o.raan, o.argp, o.nu, tc.Mu);
            rEnd = propagate(r0, v0, dt, tc.Mu);
            rExpected = kep2car(o.a, o.e, o.i, o.raan, o.argp, nu2, tc.Mu);
            tc.verifyLessThan(norm(rEnd - rExpected), 1e-3);       % [km]
        end

        function equatorialOrbitIsRejected(tc)
            [r, v] = kep2car(10000, 0.1, 0, 0, 0, 1, tc.Mu);
            tc.verifyError(@() car2kep(r, v, tc.Mu), 'car2kep:singular');
        end
    end
end

function o = randomOrbit()
o.a = 7000 + 50000*rand;
o.e = 0.01 + 0.89*rand;
o.i = deg2rad(5 + 170*rand);
o.raan = 2*pi*rand;
o.argp = 2*pi*rand;
o.nu = 2*pi*rand;
end

function d = angleDiff(x, y)
d = abs(mod(x - y + pi, 2*pi) - pi);
end

function r = propagate(r0, v0, dt, mu)
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-9);
[~, Y] = ode113(@(~, y) [y(4:6); -mu*y(1:3)/norm(y(1:3))^3], [0 dt], [r0; v0], opts);
r = Y(end, 1:3).';
end
