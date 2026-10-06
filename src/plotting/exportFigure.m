function exportFigure(fig, name, cfg)
%EXPORTFIGURE Save a figure as PNG in cfg.plot.folder when cfg.plot.export is true.

if ~cfg.plot.export
    return
end
if ~exist(cfg.plot.folder, 'dir')
    mkdir(cfg.plot.folder);
end
drawnow;
exportgraphics(fig, fullfile(cfg.plot.folder, [name '.png']), 'Resolution', 150);
end
