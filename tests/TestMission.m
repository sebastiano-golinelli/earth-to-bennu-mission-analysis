classdef TestMission < matlab.unittest.TestCase
    %TESTMISSION End-to-end checks of the Earth -> Bennu mission.
    %   Regression against reference values, plus physical checks that do
    %   not rely on the code under test: numerical propagation of the
    %   heliocentric transfer and the asymptote of the escape hyperbola.

    properties
        Cfg
        Res
    end

    methods (TestClassSetup)
        function runMissionOnce(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'config')));
            cfg = missionConfig();
            cfg.plot.enable = false;
            cfg.verbose = false;
            tc.Cfg = cfg;
            tc.Res = runMission(cfg);
        end
    end

    methods (Test)
        function regressionValues(tc)
            r = tc.Res;
            tol = {'RelTol', 1e-6};
            tc.verifyEqual(r.sc1.best.name, 'SPW pa');
            tc.verifyEqual(rad2deg(r.sc1.initial.raan), 192, 'AbsTol', 1e-9);
            tc.verifyEqual(r.sc1.best.dvTotal, 1.754656964, tol{:});
            tc.verifyEqual(r.sc1.best.time/3600, 49.14629854, tol{:});
            tc.verifyEqual(r.sc2.grid.dv, 5.126739395, tol{:});
            tc.verifyEqual(r.sc2.best.dv, 5.086991961, tol{:});
            tc.verifyEqual(r.sc2.best.tofDays, 116.293781, tol{:});
            tc.verifyEqual(r.sc3.soiTarget, 2.839145785, tol{:});
            tc.verifyEqual(r.sc3.escape.dv, 0.6824105238, tol{:});
            tc.verifyEqual(r.sc3.escape3D.dv, 0.6825460182, tol{:});
            tc.verifyEqual(r.sc3.capture.dv, 3.222851876, tol{:});
            tc.verifyEqual(r.mission.dvNonCoplanar, 5.660054857, tol{:});
        end

        function optimizersAgreeOnTheOptimum(tc)
            % Deterministic searches must find the optimum. ga is stochastic:
            % on another platform round-off can send a run to a different
            % basin, so only the validity of its result is required.
            r2 = tc.Res.sc2;
            tc.verifyEqual(r2.fmincon.dv, r2.best.dv, 'AbsTol', 1e-6);
            tc.verifyEqual(r2.multiStart.dv, r2.best.dv, 'AbsTol', 1e-6);
            tc.verifyTrue(r2.gaRefined.valid);
            tc.verifyLessThanOrEqual(r2.best.dv, r2.grid.dv);
        end

        function transferReachesTheTarget(tc)
            % Integrate the two-body problem from departure for the time of flight
            b = tc.Res.sc2.best;
            mu = tc.Cfg.const.muSun;
            opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-6);
            [~, Y] = ode113(@(~, y) [y(4:6); -mu*y(1:3)/norm(y(1:3))^3], ...
                [0 b.tofDays*86400], [b.r1; b.vT1], opts);
            tc.verifyLessThan(norm(Y(end, 1:3).' - b.r2), 1);          % [km]
            tc.verifyLessThan(norm(Y(end, 4:6).' - b.vT2), 1e-6);      % [km/s]
        end

        function escapeAsymptoteMatchesDepartureVelocity(tc)
            % Excess velocity of the escape hyperbola computed from the
            % post-burn state only: magnitude from the energy, direction
            % from the eccentricity vector.
            r3 = tc.Res.sc3;
            mu = tc.Cfg.const.muEarth;
            r = r3.escape3D.rBurn;
            v = r3.escape3D.vBurn;
            h = cross(r, v);
            eVec = cross(v, h)/mu - r/norm(r);
            e = norm(eVec);
            p = eVec/e;
            q = cross(h/norm(h), p);
            nuInf = acos(-1/e);
            vInfDir = cos(nuInf)*p + sin(nuInf)*q;
            vInfTarget = r3.vInfDepEquatorial;
            tc.verifyEqual(sqrt(norm(v)^2 - 2*mu/norm(r)), norm(vInfTarget), 'RelTol', 1e-9);
            tc.verifyLessThan(acos(min(1, dot(vInfDir, vInfTarget/norm(vInfTarget)))), 1e-6);
        end

        function designedParkingMakesEscapeCoplanar(tc)
            % With the parking orbit aligned to the asymptote, the general
            % 3D escape must cost the same as the coplanar estimate.
            r3 = tc.Res.sc3;
            tc.verifyEqual(r3.escape3D.dv, r3.escape.dv, 'AbsTol', 1e-3);   % within 1 m/s
        end

        function captureRadiiInsideBennuSoi(tc)
            T = tc.Res.sc3.captureTable;
            tc.verifyTrue(all(T.Radius_km < tc.Res.sc3.soiTarget));
            tc.verifyTrue(all(T.Radius_km > tc.Cfg.sc2.target.radius));
        end
    end
end
