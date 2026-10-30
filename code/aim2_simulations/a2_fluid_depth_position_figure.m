%A2_FLUID_DEPTH_POSITION_FIGURE  Plot peak change against fluid volume by position, depth and electrode layout.

% Aim served: Aim 2.
%
% Inputs : results/aim2/a2_sweep_index.csv (main case, 20 mm region, Maxwell Garnett mixing)
% Outputs: figures/aim2/a2_peak_by_position_depth.png, figures/aim2/a2_peak_by_layout.png
% Imogen Ryan - Honours thesis.

clear; close all force;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
T = readtable(fullfile(root,'results','aim2','a2_sweep_index.csv'), 'TextType','string');
T = T(T.scenario == "B" & T.rule == "maxwell_garnett", :);     % main case, 558 runs
T.peak_uV = T.peak_dV_mV * 1e3;                                 % peak change (uV)
T.pct     = T.burden_frac * 100;                               % fluid in region (% v/v)
levels    = unique(T.pct);                                     % 0.25 0.5 1 2 5 10
vols      = levels/100*pi/6*20^3;                              % extra fluid (mm^3)
T.vol     = T.burden_frac*pi/6*20^3;                           % extra fluid per run (mm^3)
tick_v    = [0 vols(4:6).'];                                     % 0, 2%, 5%, 10% (smaller levels too close to label)
tick_lab  = [{'0\newline(0%)'}, compose('%.1f\\newline(%g%%)', vols(4:6), levels(4:6)).'];
xfit      = [0 450];                                           % x range of the best-fit lines (mm^3)

% Graph 1: one panel per scalp point, one line per depth (19 electrodes)
M = T(T.montage == "m19_1020", :);
points = ["Cz","Fz","C3","Pz","T7","O1"];
where  = ["top of head","front, middle","left, top","back, top","left side, low","back, low"];
lit    = ["above the centrum semiovale (Alzheimer's)", "above frontal white matter (brain injury)", "", "", "", ""];
depths = 25:5:55;
cmap   = flipud(parula(numel(depths)+2)); cmap = cmap(2:end-1,:);  % light = shallow

f = figure('Units','centimeters','Position',[1 1 17 17],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
tl = tiledlayout(f, 3, 2, 'Padding','compact', 'TileSpacing','compact');   % three rows of two so the labels fit at 10 pt
for p = 1:numel(points)
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'off'); box(ax,'off')
    for k = 1:numel(depths)
        m = M.anchor == points(p) & M.depth_mm == depths(k);
        if ~any(m), continue; end
        s = M(m,:);
        slope = sum(s.vol.*s.peak_uV) / sum(s.vol.^2);         % least-squares line through zero
        plot(ax, xfit, slope*xfit, '--', 'Color', cmap(k,:), 'LineWidth', 1, 'HandleVisibility','off');
        plot(ax, s.vol, s.peak_uV, 'o', 'Color', cmap(k,:), ...
            'MarkerFaceColor', cmap(k,:), 'MarkerSize', 4, ...
            'DisplayName', sprintf('%d mm deep', depths(k)));
    end
    xlim(ax,[0 450]); ylim(ax,[0 2.6]);
    xticks(ax,tick_v); xticklabels(ax,tick_lab); ax.TickLabelInterpreter = 'tex'; ax.XTickLabelRotation = 0;
    if p >= 5, xlabel(ax,'Extra fluid in a 20 mm region of brain (mm^3)'); end       % bottom row only
    if mod(p,2) == 1, ylabel(ax,{'Peak change (\muV)','with fluid − healthy'}); end   % left column only
    title(ax, sprintf('%s (%s)', points(p), where(p)));
    if lit(p) ~= "", ax.Title.Color = [0.75 0.10 0.10]; ax.LineWidth = 1.5; end   % literature-review positions
end
% one legend for all panels, listing every depth
hd = gobjects(numel(depths),1);
for k = 1:numel(depths)
    hd(k) = plot(ax, NaN, NaN, '--o', 'Color', cmap(k,:), 'MarkerFaceColor', cmap(k,:), ...
        'LineWidth', 1, 'DisplayName', sprintf('%d mm', depths(k)));
end
lg = legend(hd, 'Orientation','horizontal'); lg.Layout.Tile = 'south';
title(lg, 'Depth of the fluid region below the scalp point'); legend(ax,'boxoff');
exportgraphics(f, fullfile(root,'figures','aim2','a2_peak_by_position_depth.png'), 'Resolution', 300);

% Graph 2: the three electrode layouts, median and SD over all positions
lay  = ["m19_1020","m32_1010","m16_twoplane"];
name = ["19 electrodes (10-20)","32 electrodes (10-10)","16 electrodes (two planes)"];
col  = [0.20 0.40 0.70; 0.93 0.60 0.20; 0.30 0.65 0.45];
g = figure('Units','centimeters','Position',[1 1 14 10],'Color','w','Visible','off'); theme(g,'light');
set(g, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
ax = axes(g); hold(ax,'on'); grid(ax,'off'); box(ax,'off')
for q = 1:numel(lay)
    L = T(T.montage == lay(q), :);
    b = unique(L.burden_frac);
    v = b*pi/6*20^3 + (q-2)*4;                            % fluid (mm^3), nudged 4 mm^3 so the bars do not overlap
    med = zeros(numel(b),1); sd = med;
    for j = 1:numel(b)
        y = L.peak_uV(L.burden_frac == b(j));
        med(j) = median(y); sd(j) = std(y);
    end
    vb = b*pi/6*20^3;                                     % fluid without the nudge (mm^3)
    slope = sum(vb.*med) / sum(vb.^2);                    % least-squares line through zero, fitted to the medians
    plot(ax, xfit, slope*xfit, '--', 'Color', col(q,:), 'LineWidth', 1.2, 'HandleVisibility','off');
    errorbar(ax, v, med, sd, 'o', 'Color', col(q,:), 'MarkerFaceColor', col(q,:), ...
        'LineWidth', 1.4, 'CapSize', 5, 'DisplayName', name(q));
end
xticks(ax,tick_v); xticklabels(ax,tick_lab); ax.TickLabelInterpreter = 'tex'; xlim(ax,[0 450]); ylim(ax,[0 2.7]);
xlabel(ax,'Extra fluid in a 20 mm region of brain (mm^3)'); ylabel(ax,'Peak change, with fluid − healthy (\muV)');
lg2 = legend(ax,'Location','northwest'); legend(ax,'boxoff');
exportgraphics(g, fullfile(root,'figures','aim2','a2_peak_by_layout.png'), 'Resolution', 300);
