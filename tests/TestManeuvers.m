classdef TestManeuvers < matlab.unittest.TestCase
    %TESTMANEUVERS Geometric checks of the impulsive maneuvers: the impulse is
    %   applied at a point that belongs to both the old and the new orbit,
    %   and its magnitude equals the velocity difference there.

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
        function bitangentMatchesHohmann(tc)
            r1 = 7000; r2 = 42164;
            [dv1, dv2, dt] = bitangentTransfer(r1, 0, r2, 0, 'pa', 0, tc.Mu);
            aT = (r1 + r2)/2;
            dvHohmann = sqrt(tc.Mu/r1)*(sqrt(2*r2/(r1 + r2)) - 1) + sqrt(tc.Mu/r2)*(1 - sqrt(2*r1/(r1 + r2)));
            tc.verifyEqual(dv1 + dv2, dvHohmann, 'RelTol', 1e-12);
            tc.verifyEqual(dt, pi*sqrt(aT^3/tc.Mu), 'RelTol', 1e-12);
        end

        function planeChangeKeepsPositionAndRotatesPlane(tc)
            % all four sign combinations of the RAAN and inclination changes
            rng(3);
            for k = 1:200
                a = 10000 + 30000*rand; e = 0.8*rand;
                iI = deg2rad(5 + 80*rand); iF = deg2rad(5 + 80*rand);
                raanI = 2*pi*rand; raanF = 2*pi*rand; argpI = 2*pi*rand;
                [dv, argpF, nu] = planeChange(a, e, iI, raanI, argpI, iF, raanF, tc.Mu);
                [rOld, vOld] = kep2car(a, e, iI, raanI, argpI, nu, tc.Mu);
                [rNew, vNew] = kep2car(a, e, iF, raanF, argpF, nu, tc.Mu);
                tc.verifyLessThan(norm(rNew - rOld), 1e-6*norm(rOld));
                tc.verifyEqual(norm(vNew - vOld), dv, 'RelTol', 1e-9);
                tc.verifyLessThanOrEqual(cos(nu), 1e-12);           % node on the apoapsis side
            end
        end

        function periapsisChangeKeepsPosition(tc)
            rng(4);
            for k = 1:100
                a = 10000 + 30000*rand; e = 0.05 + 0.85*rand;
                i = deg2rad(5 + 80*rand); raan = 2*pi*rand;
                argpI = 2*pi*rand; argpF = 2*pi*rand;
                [dv, nuI, nuF] = periapsisChange(a, e, argpI, argpF, tc.Mu);
                for j = 1:2
                    [rOld, vOld] = kep2car(a, e, i, raan, argpI, nuI(j), tc.Mu);
                    [rNew, vNew] = kep2car(a, e, i, raan, argpF, nuF(j), tc.Mu);
                    tc.verifyLessThan(norm(rNew - rOld), 1e-6*norm(rOld));
                    tc.verifyEqual(norm(vNew - vOld), dv, 'RelTol', 1e-9);
                end
            end
        end

        function bitangentWithPlaneChangeIsContinuous(tc)
            oI = struct('a', 24400, 'e', 0.73, 'i', deg2rad(27), 'raan', 0.4, 'argp', 3.1, 'nu', 0);
            oF = struct('a', 36700, 'e', 0.81, 'i', deg2rad(35), 'raan', 1.0, 'argp', 0, 'nu', 0);
            nTested = 0;
            for type = {'pa', 'pp', 'aa'}
                m = bitangentWithPlaneChange(oI, oF, type{1}, tc.Mu);
                if ~m.feasible
                    continue
                end
                nTested = nTested + 1;
                % the plane change happens at a point common to both transfer planes
                [rPre, vPre] = kep2car(m.aT, m.eT, oI.i, oI.raan, m.argpBefore, m.nuNode, tc.Mu);
                [rPost, vPost] = kep2car(m.aT, m.eT, oF.i, oF.raan, m.argpAfter, m.nuNode, tc.Mu);
                tc.verifyLessThan(norm(rPost - rPre), 1e-6*norm(rPre));
                tc.verifyEqual(norm(vPost - vPre), m.dvPlane, 'RelTol', 1e-9);
                % the transfer starts on the initial orbit and ends on the final one
                rDep = kep2car(oI.a, oI.e, oI.i, oI.raan, oI.argp, apsisAnomaly(type{1}(1)), tc.Mu);
                rStart = kep2car(m.aT, m.eT, oI.i, oI.raan, m.argpBefore, 0, tc.Mu);
                tc.verifyLessThan(norm(rStart - rDep), 1e-6*norm(rDep));
                rArr = kep2car(oF.a, oF.e, oF.i, oF.raan, m.argpFinal, apsisAnomaly(type{1}(2)), tc.Mu);
                rEnd = kep2car(m.aT, m.eT, oF.i, oF.raan, m.argpAfter, pi, tc.Mu);
                tc.verifyLessThan(norm(rEnd - rArr), 1e-6*norm(rArr));
            end
            tc.verifyGreaterThan(nTested, 0);
        end
    end
end
