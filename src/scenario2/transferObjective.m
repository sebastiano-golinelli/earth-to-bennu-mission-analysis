function J = transferObjective(x, cfg)
%TRANSFEROBJECTIVE Objective for the optimizers: cost, or a large penalty if invalid.

sol = evaluateTransfer(x, cfg);
if sol.valid
    J = sol.J;
else
    J = 1e9 + sol.penalty;
end
if ~isfinite(J)
    J = 1e12;
end
end
