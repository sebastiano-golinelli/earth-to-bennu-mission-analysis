function plotScenario3(cfg, r3)
%PLOTSCENARIO3 Figures of the Earth escape and of the capture at Bennu.

c = cfg.const;
park = r3.parking;
H = r3.escape3D;

% 3D escape from the parking orbit
fig = figure('Name', 'SC3 - Escape hyperbola', 'Color', 'w');
hold on; grid on; axis equal; box on;
plotEarth(c.rEarth);
plotOrbitArc(park.a, park.e, park.i, park.raan, park.argp, 0, 2*pi, 2e-3, [0.85 0.33 0.10], 'Parking orbit', '-', 1.5);
rMax = 2.5*park.ra;
nuH = linspace(H.nuBurn, 0.98*H.nuInf, 2000);
P = conicPoints(H.a, H.e, H.i, H.raan, H.argp, nuH);
P = P(:, vecnorm(P) <= rMax);
plot3(P(1, :), P(2, :), P(3, :), 'r', 'LineWidth', 2, 'DisplayName', 'Escape hyperbola');
plot3(H.rBurn(1), H.rBurn(2), H.rBurn(3), 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 12, 'DisplayName', 'Escape burn');
d = r3.vInfDepEquatorial/norm(r3.vInfDepEquatorial)*rMax;
quiver3(0, 0, 0, d(1), d(2), d(3), 0, 'm', 'LineWidth', 1.6, 'MaxHeadSize', 0.3, 'DisplayName', 'Departure v_\infty direction');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title(sprintf('Escape: \\Deltav = %.4f km/s, v_\\infty = %.3f km/s', H.dv, r3.vInfDep));
legend('Location', 'bestoutside'); view(35, 25);
exportFigure(fig, 'sc3_escape', cfg);

% Escape cost versus burn point along the parking orbit
fig = figure('Name', 'SC3 - Burn point', 'Color', 'w');
plot(rad2deg(r3.burnSweep.nu), r3.burnSweep.dv, 'LineWidth', 1.6, 'DisplayName', '3D escape from this point'); hold on; grid on;
yline(r3.escape.dv, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Coplanar estimate (tangential burn at periapsis)');
plot(rad2deg(H.nuBurnOnParking), H.dv, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 12, 'DisplayName', 'Selected burn point');
xlim([0 360]); xticks(0:60:360);
set(gca, 'YScale', 'log');
xlabel('True anomaly of the burn point on the parking orbit [deg]');
ylabel('Escape \Deltav [km/s]');
title({'Non-coplanar escape: cost of the burn point', ...
    '(gaps: the escape hyperbola would pass below the Earth''s surface)'});
legend('Location', 'northwest');
exportFigure(fig, 'sc3_burn_point', cfg);

% Capture at Bennu: validity of the two-body model
T = r3.captureTable;
fig = figure('Name', 'SC3 - Capture at Bennu', 'Color', 'w');
t = tiledlayout(fig, 'flow', 'TileSpacing', 'compact');
nexttile(t); hold on; grid on; box on;
plot(T.Radius_km, T.RadiusOverSOI, 'o-', 'LineWidth', 1.4, 'DisplayName', 'r / r_{SOI}');
plot(T.Radius_km, T.SRPRatio, 's-', 'LineWidth', 1.4, 'DisplayName', 'Solar radiation pressure / gravity');
plot(T.Radius_km, T.TideRatio, 'd-', 'LineWidth', 1.4, 'DisplayName', 'Solar tide / gravity');
yline(0.1, 'k:', '10%', 'HandleVisibility', 'off');
set(gca, 'YScale', 'log');
xlabel('Capture orbit radius [km]'); ylabel('Ratio [-]');
title('Is the two-body capture model valid?');
legend('Location', 'southeast');
nexttile(t); hold on; grid on; box on;
plot(T.Radius_km, 1000*(r3.vInfArrMag - T.DeltaVCircular_kms), 'o-', 'LineWidth', 1.4);
xlabel('Capture orbit radius [km]'); ylabel('v_\infty - \Deltav_{capture} [m/s]');
title('Saving due to Bennu''s gravity');
title(t, sprintf('Capture at Bennu (SOI radius %.2f km, v_\\infty = %.3f km/s)', r3.soiTarget, r3.vInfArrMag));
exportFigure(fig, 'sc3_capture', cfg);
end
