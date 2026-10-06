function h = plotOrbitArc(a, e, i, raan, argp, nu0, nu1, dnu, color, name, style, width)
%PLOTORBITARC Plot an orbit arc from true anomaly nu0 to nu1 on the current axes.
%   h = PLOTORBITARC(a, e, i, raan, argp, nu0, nu1, dnu, color, name, style, width)
%   draws the arc with step dnu [rad], a circle at the start and a square
%   at the end. name is used in the legend.

if nargin < 11 || isempty(style), style = '-'; end
if nargin < 12 || isempty(width), width = 1.5; end
nu = nu0:dnu:nu1;
if isempty(nu) || nu(end) < nu1, nu = [nu, nu1]; end
P = conicPoints(a, e, i, raan, argp, nu);
h = plot3(P(1, :), P(2, :), P(3, :), style, 'Color', color, 'LineWidth', width, 'DisplayName', name);
plot3(P(1, 1), P(2, 1), P(3, 1), 'o', 'Color', color, 'MarkerFaceColor', color, ...
    'MarkerSize', 4, 'HandleVisibility', 'off');
plot3(P(1, end), P(2, end), P(3, end), 's', 'Color', color, 'MarkerFaceColor', color, ...
    'MarkerSize', 4, 'HandleVisibility', 'off');
end
