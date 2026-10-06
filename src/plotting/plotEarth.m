function plotEarth(radius)
%PLOTEARTH Earth with topographic texture centred at the origin of the current axes.
%   PLOTEARTH(radius) uses MATLAB's 'topo' dataset (longitude 0..360 deg).
%   The globe is slightly transparent, so that orbit points behind it stay visible.

if nargin < 1, radius = 6378.137; end
S = load('topo', 'topo');
[nLat, nLon] = size(S.topo);
[lon, lat] = meshgrid(linspace(0, 2*pi, nLon), linspace(-pi/2, pi/2, nLat));
surf(radius*cos(lat).*cos(lon), radius*cos(lat).*sin(lon), radius*sin(lat), S.topo, ...
    'FaceColor', 'texturemap', 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');
colormap(gca, earthColormap());
clim([-8000 8000]);
end

function cmap = earthColormap()
sea      = [linspace(0.00, 0.10, 90)' linspace(0.20, 0.50, 90)' linspace(0.50, 1.00, 90)'];
land     = [linspace(0.10, 0.25, 70)' linspace(0.45, 0.80, 70)' linspace(0.10, 0.20, 70)'];
mountain = [linspace(0.35, 0.60, 50)' linspace(0.25, 0.45, 50)' linspace(0.15, 0.25, 50)'];
ice      = [linspace(0.85, 1.00, 46)' linspace(0.85, 1.00, 46)' linspace(0.85, 1.00, 46)'];
cmap = [sea; land; mountain; ice];
end
