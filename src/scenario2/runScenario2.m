function r2 = runScenario2(cfg)
%RUNSCENARIO2 Direct two-impulse heliocentric transfer from Earth to the target.
%   r2 = RUNSCENARIO2(cfg) optimizes x = [nuDep, nuArr, argpT] with
%     1) a grid search,
%     2) fmincon (SQP) started from the best grid point,
%     3) ga, several seeds, best run refined by fmincon,
%     4) MultiStart (many fmincon runs from random points),
%   and returns every result plus the best valid one in r2.best.
%   Steps 3 and 4 need the Global Optimization Toolbox and are skipped
%   without it.

s = cfg.sc2;
obj = @(x) transferObjective(x, cfg);
hasGlobalOpt = exist('ga', 'file') == 2 && ~isempty(which('MultiStart'));

% ---- 1) Grid search -------------------------------------------------------
nuDepGrid = gridVector(s.grid(1));
nuArrGrid = gridVector(s.grid(2));
argpGrid  = gridVector(s.grid(3));
n = prod(s.grid);
rec = struct('x', nan(n, 3), 'dv', nan(n, 1), 'tofDays', nan(n, 1), 'arcDeg', nan(n, 1), ...
    'valid', false(n, 1));
best = evaluateTransfer([0 0 0], cfg);
best.J = Inf;
k = 0;
for i1 = 1:s.grid(1)
    for i2 = 1:s.grid(2)
        for i3 = 1:s.grid(3)
            k = k + 1;
            x = [nuDepGrid(i1), nuArrGrid(i2), argpGrid(i3)];
            sol = evaluateTransfer(x, cfg);
            rec.x(k, :) = x;
            rec.dv(k) = sol.dv;
            rec.tofDays(k) = sol.tofDays;
            rec.arcDeg(k) = sol.arcDeg;
            rec.valid(k) = sol.valid;
            if sol.valid && sol.J < best.J
                best = sol;
            end
        end
    end
end
if ~any(rec.valid)
    error('runScenario2:noSolution', 'No valid transfer on the grid.');
end
best.method = "Grid search";
r2.grid = best;
r2.gridRecords = rec;
r2.gridVectors = {nuDepGrid, nuArrGrid, argpGrid};

% ---- 2) fmincon from the best grid point ----------------------------------
% No bounds: the objective is periodic in all three angles.
histF = struct('iter', [], 'J', []);
optF = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxIterations', s.fmincon.maxIterations, 'MaxFunctionEvaluations', 1e4, ...
    'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-12, 'OutputFcn', @recordFmincon);
xF = fmincon(obj, r2.grid.x, [], [], [], [], [], [], [], optF);
r2.fmincon = finalize(xF, "fmincon");
r2.fminconHistory = histF;

% Same local refinement, without the iteration history (used after ga)
optPolish = optimoptions(optF, 'OutputFcn', []);

% ---- 3) ga, several seeds -------------------------------------------------
r2.ga = finalize([NaN NaN NaN], "GA");
r2.gaRefined = finalize([NaN NaN NaN], "GA + fmincon");
r2.gaRuns = [];
r2.gaHistory = struct('generation', [], 'J', []);
if hasGlobalOpt
    optGA = optimoptions('ga', 'Display', 'off', 'MaxGenerations', s.ga.generations, ...
        'PopulationSize', s.ga.population, 'FunctionTolerance', 1e-8, 'OutputFcn', @recordGa);
    bestJ = Inf;
    for seed = s.ga.seeds
        rng(seed);
        histG = struct('generation', [], 'J', []);
        [xG, JG] = ga(obj, 3, [], [], [], [], zeros(1, 3), 2*pi*ones(1, 3), [], optGA);
        r2.gaRuns = [r2.gaRuns; seed, JG];
        if JG < bestJ
            bestJ = JG;
            r2.ga = finalize(xG, "GA");
            r2.gaHistory = histG;
        end
    end
    r2.gaRefined = finalize(fmincon(obj, r2.ga.x, [], [], [], [], [], [], [], optPolish), "GA + fmincon");
end

% ---- 4) MultiStart ----------------------------------------------------------
r2.multiStart = finalize([NaN NaN NaN], "MultiStart");
r2.localMinima = table();
if hasGlobalOpt && s.multiStart.points > 0
    rng(s.multiStart.seed);
    % bounds only define where the random start points are drawn
    prob = createOptimProblem('fmincon', 'objective', obj, 'x0', r2.grid.x, ...
        'lb', -pi*ones(1, 3), 'ub', 3*pi*ones(1, 3), 'options', optPolish);
    [xM, ~, ~, ~, sols] = run(MultiStart('Display', 'off'), prob, s.multiStart.points);
    r2.multiStart = finalize(xM, "MultiStart");
    r2.localMinima = localMinimaTable(sols, cfg);
end

% ---- Best valid solution ----------------------------------------------------
cand = [r2.grid, r2.fmincon, r2.gaRefined, r2.multiStart];
cand = cand([cand.valid]);
[~, kBest] = min([cand.J]);
r2.best = cand(kBest);
m = [r2.grid, r2.fmincon, r2.ga, r2.gaRefined, r2.multiStart];
r2.comparison = table([m.method].', [m.dv].', [m.dv1].', [m.dv2].', [m.tofDays].', ...
    [m.arcDeg].', [m.valid].', 'VariableNames', ...
    {'Method', 'DeltaV_kms', 'DeltaV1_kms', 'DeltaV2_kms', 'TOF_days', 'Arc_deg', 'Valid'});

    function sol = finalize(x, method)
        sol = evaluateTransfer(x, cfg);
        sol.method = method;
    end

    function stop = recordFmincon(~, optimValues, state)
        stop = false;
        if strcmp(state, 'iter')
            histF.iter(end+1, 1) = optimValues.iteration;
            histF.J(end+1, 1) = optimValues.fval;
        end
    end

    function [state, options, changed] = recordGa(options, state, flag)
        changed = false;
        if strcmp(flag, 'iter') || strcmp(flag, 'init')
            histG.generation(end+1, 1) = state.Generation;
            histG.J(end+1, 1) = min(state.Score);
        end
    end
end

function v = gridVector(n)
v = linspace(0, 2*pi, n + 1);
v(end) = [];
end

function T = localMinimaTable(sols, cfg)
% Distinct local minima found by MultiStart, with how many starts reached each one
nSol = numel(sols);
dv = nan(nSol, 1); tof = dv; arc = dv; starts = dv; valid = false(nSol, 1);
for k = 1:nSol
    s = evaluateTransfer(sols(k).X, cfg);
    dv(k) = s.dv; tof(k) = s.tofDays; arc(k) = s.arcDeg; valid(k) = s.valid;
    starts(k) = numel(sols(k).X0);
end
T = table(dv, tof, arc, starts, 'VariableNames', {'DeltaV_kms', 'TOF_days', 'Arc_deg', 'Starts'});
T = T(valid, :);
% The same transfer can be found with angles shifted by 2*pi: merge duplicates
[~, first, group] = unique([round(T.DeltaV_kms, 5), round(T.TOF_days, 2)], 'rows', 'stable');
nStarts = accumarray(group, T.Starts);
T = T(first, :);
T.Starts = nStarts;
T = sortrows(T, 'DeltaV_kms');
end
