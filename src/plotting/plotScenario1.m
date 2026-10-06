function plotScenario1(cfg, r1)
%PLOTSCENARIO1 Figures of the geocentric transfer.

dnu = 2e-3;
oI = r1.initial;
oF = r1.final;

% Launch and parking orbits
fig = figure('Name', 'SC1 - Launch and parking orbits', 'Color', 'w');
hold on; grid on; axis equal;
plotEarth(cfg.const.rEarth);
plotOrbitArc(oI.a, oI.e, oI.i, oI.raan, oI.argp, 0, 2*pi, dnu, [0 0.45 0.74], 'Launch orbit (GTO)', '-', 1.6);
plotOrbitArc(oF.a, oF.e, oF.i, oF.raan, oF.argp, 0, 2*pi, dnu, [0.85 0.33 0.10], 'Parking orbit', '-', 1.6);
r0 = r1.initialState.r;
plot3(r0(1), r0(2), r0(3), 'kd', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'DisplayName', 'Spacecraft at start');
P = conicPoints(oF.a, oF.e, oF.i, oF.raan, oF.argp, 0);
plot3(P(1), P(2), P(3), 'kp', 'MarkerFaceColor', 'r', 'MarkerSize', 11, 'DisplayName', 'Escape burn point');
labelAxes3D('Launch orbit and parking orbit (equatorial frame)');
exportFigure(fig, 'sc1_orbits', cfg);

% Launch window: geocentric cost as a function of the launch RAAN
if ~isempty(r1.launchSweep)
    T = r1.launchSweep;
    fig = figure('Name', 'SC1 - Launch window', 'Color', 'w');
    plot(T.LaunchRAAN_deg, T.DeltaV_kms, 'LineWidth', 1.6); hold on; grid on;
    plot(rad2deg(oI.raan), r1.best.dvTotal, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
    xlim([0 360]); xticks(0:60:360);
    xlabel('RAAN of the launch orbit [deg] (set by the launch time)');
    ylabel('Geocentric \Deltav [km/s]');
    title('Launch window: cheapest strategy for each launch RAAN');
    legend('Best strategy', sprintf('Selected: %.0f deg, %.3f km/s', rad2deg(oI.raan), r1.best.dvTotal), ...
        'Location', 'north');
    exportFigure(fig, 'sc1_launch_window', cfg);
end

% Strategy trade-off: cost and duration of every strategy (rows) and
% bitangent type (columns); cases are ordered strategy first, then type
nS = numel(cfg.sc1.strategies);
nT = numel(cfg.sc1.types);
dv = reshape(r1.table.DeltaV_kms, nT, nS).';
time = reshape(r1.table.Time_h, nT, nS).';
fig = figure('Name', 'SC1 - Strategy trade-off', 'Color', 'w');
h = heatmap(fig, cfg.sc1.types, cfg.sc1.strategies, dv, 'Position', [0.08 0.12 0.36 0.72], ...
    'CellLabelFormat', '%.3f', 'ColorbarVisible', 'off');
h.Title = [char(916) 'v [km/s]'];
h.XLabel = 'Bitangent (departure, arrival apsis)';
h.YLabel = 'Strategy';
h = heatmap(fig, cfg.sc1.types, cfg.sc1.strategies, time, 'Position', [0.58 0.12 0.36 0.72], ...
    'CellLabelFormat', '%.1f', 'ColorbarVisible', 'off');
h.Title = 'Transfer time [h]';
h.XLabel = 'Bitangent (departure, arrival apsis)';
annotation(fig, 'textbox', [0 0.9 1 0.1], 'String', sprintf( ...
    'Geocentric strategies (S = shape, P = plane, W = periapsis) - selected: %s', r1.best.name), ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 11);
exportFigure(fig, 'sc1_strategies', cfg);

% Maneuver sequence of the selected strategy (optionally all of them)
if cfg.plot.allSC1Cases
    idx = find([r1.cases.feasible]);
else
    idx = r1.selected;
end
for k = idx
    c = r1.cases(k);
    fig = figure('Name', ['SC1 - ' c.name], 'Color', 'w');
    hold on; grid on; axis equal;
    plotEarth(cfg.const.rEarth);
    plotStrategySequence(c, oI, oF, dnu);
    labelAxes3D(sprintf('Strategy %s: \\Deltav = %.3f km/s, %.1f h', c.name, c.dvTotal, min(c.dtTotal)/3600));
    if k == r1.selected
        exportFigure(fig, 'sc1_selected_strategy', cfg);
    end
end
end

function labelAxes3D(ttl)
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title(ttl);
legend('Location', 'bestoutside');
view(3);
end
