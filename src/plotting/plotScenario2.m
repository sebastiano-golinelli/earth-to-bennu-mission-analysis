function plotScenario2(cfg, r2)
%PLOTSCENARIO2 Figures of the heliocentric transfer and of the optimization.

AU = cfg.const.AU;
rec = r2.gridRecords;
best = r2.best;

% Grid search: minimum dv over argpT for each departure/arrival pair
[nuDep, nuArr, ~] = r2.gridVectors{:};
n = cfg.sc2.grid;
dv = rec.dv;
dv(~rec.valid) = NaN;
dvMin = squeeze(min(reshape(dv, n(3), n(2), n(1)), [], 1)).';    % n(1) x n(2)
fig = figure('Name', 'SC2 - Grid search', 'Color', 'w');
imagesc(rad2deg(nuArr), rad2deg(nuDep), dvMin, 'AlphaData', double(~isnan(dvMin)));
set(gca, 'YDir', 'normal'); hold on; colorbar;
clim([min(dvMin(:)), min(dvMin(:)) + 10]);
plot(rad2deg(best.x(2)), rad2deg(best.x(1)), 'wp', 'MarkerFaceColor', 'r', 'MarkerSize', 14);
xlabel('Arrival true anomaly on Bennu''s orbit [deg]');
ylabel('Departure true anomaly on Earth''s orbit [deg]');
title('Grid search: minimum \Deltav over \omega_T [km/s]');
exportFigure(fig, 'sc2_grid_search', cfg);

% Cost versus time of flight: grid points, local minima, best solution
fig = figure('Name', 'SC2 - Delta-v vs time of flight', 'Color', 'w');
hold on; grid on; box on;
v = rec.valid;
scatter(rec.tofDays(v), rec.dv(v), 5, rec.arcDeg(v), 'filled', 'DisplayName', 'Grid points');
cb = colorbar; cb.Label.String = 'Transfer angle [deg]';
if ~isempty(r2.localMinima)
    L = r2.localMinima;
    plot(L.TOF_days, L.DeltaV_kms, 'ko', 'MarkerSize', 7, 'LineWidth', 1.1, 'DisplayName', 'MultiStart local minima');
end
plot(best.tofDays, best.dv, 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 16, 'DisplayName', 'Best transfer');
ylim([0.95*best.dv, best.dv + 5]);
xlabel('Time of flight [days]'); ylabel('\Deltav_1 + \Deltav_2 [km/s]');
title('Heliocentric transfers: cost versus time of flight');
legend('Location', 'northeast');
exportFigure(fig, 'sc2_pareto', cfg);

% Best transfer: top view and 3D view
nu = linspace(0, 2*pi, 600);
oE = cfg.sc2.earth; oT = cfg.sc2.target;
PE = conicPoints(oE.a, oE.e, oE.i, oE.raan, oE.argp, nu)/AU;
PB = conicPoints(oT.a, oT.e, oT.i, oT.raan, oT.argp, nu)/AU;
PT = conicPoints(best.aT, best.eT, best.iT, best.raanT, best.argpT, ...
    best.nu1T + linspace(0, deg2rad(best.arcDeg), 400))/AU;
r1 = best.r1/AU; r2pos = best.r2/AU;
cE = [0 0.45 0.74]; cB = [0.47 0.47 0.47];
fig = figure('Name', 'SC2 - Best transfer', 'Color', 'w');
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% top view of the ecliptic plane
ax = nexttile(t); hold on; grid on; box on; axis equal;
plot(PE(1, :), PE(2, :), 'Color', cE, 'LineWidth', 1.4, 'DisplayName', 'Earth orbit');
plot(PB(1, :), PB(2, :), 'Color', cB, 'LineWidth', 1.4, 'DisplayName', 'Bennu orbit');
plot(PT(1, :), PT(2, :), 'r', 'LineWidth', 2.4, 'DisplayName', 'Transfer');
plot(0, 0, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.8 0], 'MarkerEdgeColor', 'k', 'DisplayName', 'Sun');
plot(r1(1), r1(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', cE, 'MarkerEdgeColor', 'k', 'DisplayName', 'Departure');
plot(r2pos(1), r2pos(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', cB, 'MarkerEdgeColor', 'k', 'DisplayName', 'Arrival');
xlabel('x [AU]'); ylabel('y [AU]');
title('Ecliptic plane, top view');
lg = legend(ax, 'NumColumns', 6);
lg.Layout.Tile = 'south';
% height above the ecliptic versus ecliptic longitude
nexttile(t); hold on; grid on; box on;
[lon, z] = longitudeHeight(PE); plot(lon, z, 'Color', cE, 'LineWidth', 1.4);
[lon, z] = longitudeHeight(PB); plot(lon, z, 'Color', cB, 'LineWidth', 1.4);
[lon, z] = longitudeHeight(PT); plot(lon, z, 'r', 'LineWidth', 2.4);
[lon, z] = longitudeHeight(r1); plot(lon, z, 'o', 'MarkerSize', 8, 'MarkerFaceColor', cE, 'MarkerEdgeColor', 'k');
[lon, z] = longitudeHeight(r2pos); plot(lon, z, 'o', 'MarkerSize', 8, 'MarkerFaceColor', cB, 'MarkerEdgeColor', 'k');
xlim([0 360]); xticks(0:60:360);
xlabel('Ecliptic longitude [deg]'); ylabel('Height above the ecliptic [AU]');
title('Out-of-plane motion: arrival at Bennu''s node');
title(t, sprintf('Best transfer: \\Deltav = %.3f km/s, TOF = %.0f days, transfer angle %.0f deg', ...
    best.dv, best.tofDays, best.arcDeg));
exportFigure(fig, 'sc2_best_transfer', cfg);

% Optimizers
fig = figure('Name', 'SC2 - Optimizers', 'Color', 'w');
t = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t);
plot(r2.fminconHistory.iter, r2.fminconHistory.J, 'o-', 'LineWidth', 1.3); grid on;
xlabel('Iteration'); ylabel('\Deltav [km/s]'); title('fmincon from the best grid point');
nexttile(t);
if ~isempty(r2.gaHistory.generation)
    plot(r2.gaHistory.generation, r2.gaHistory.J, 'LineWidth', 1.4); grid on;
    ylim([best.dv - 0.05, min(max(r2.gaHistory.J), best.dv + 2)]);
end
xlabel('Generation'); ylabel('Best \Deltav [km/s]'); title('ga, best of the seeds');
nexttile(t); hold on; grid on; box on;
if ~isempty(r2.gaRuns)
    bar(categorical("seed " + string(r2.gaRuns(:, 1))), r2.gaRuns(:, 2));
    yline(best.dv, 'r--', 'Global optimum', 'LineWidth', 1.3, 'LabelHorizontalAlignment', 'left');
    ylim([best.dv - 0.1, max(r2.gaRuns(:, 2)) + 0.1]);
end
ylabel('Final \Deltav [km/s]'); title('ga: result of each seed');
exportFigure(fig, 'sc2_optimizers', cfg);
end

function [lon, z] = longitudeHeight(P)
% ecliptic longitude [deg] and height [same unit as P] of 3xN points,
% with NaN where the longitude wraps around 360 deg (no spurious lines)
lon = mod(rad2deg(atan2(P(2, :), P(1, :))), 360);
z = P(3, :);
jump = find(abs(diff(lon)) > 180);
for k = numel(jump):-1:1
    lon = [lon(1:jump(k)), NaN, lon(jump(k)+1:end)];
    z = [z(1:jump(k)), NaN, z(jump(k)+1:end)];
end
end
