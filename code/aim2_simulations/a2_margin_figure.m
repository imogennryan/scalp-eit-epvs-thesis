%A2_MARGIN_FIGURE  Plot the margin above the detection threshold against the extra fluid, for four machines.

% Aim served: Aim 2.
%
% Inputs : results/aim2/a2_noise_summary.csv (from a2_add_noise.m). Main case
%          (20 mm region, Maxwell Garnett), all three electrode layouts.
% Outputs: figures/aim2/a2_margin_by_fluid.png
% Imogen Ryan - Honours thesis.

clear; close all force;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
T = readtable(fullfile(root,'results','aim2','a2_noise_summary.csv'), 'TextType','string');
T = T(T.scenario == "B" & T.rule == "maxwell_garnett", :);     % main case

levels = unique(T.burden_frac)*100;                  % 0.25 0.5 1 2 5 10 (% v/v)
[~, T.lvl] = ismember(T.burden_frac*100, levels);   % evenly spaced position 1-6
machines = [45 48 77 95];
names = ["Clinical, higher noise (45 dB)","Clinical, lower noise (48 dB)","Research, higher noise (77 dB)","Research, lowest noise (95 dB)"];
col   = [0.80 0.30 0.30; 0.95 0.55 0.45; 0.35 0.55 0.85; 0.15 0.30 0.65];

f = figure('Units','centimeters','Position',[1 1 15 10],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
ax = axes(f); hold(ax,'on'); grid(ax,'on'); box(ax,'off')

% position on the evenly spaced axis for a fluid % between two levels (log interpolation)
pos = @(pct) interp1(log(levels), 1:numel(levels), log(pct));
x77 = pos(100*254.5/4188.79); x95 = pos(100*31.7/4188.79);   % crossing points (a2_results_summary)
yl = [-75 45];
patch(ax, [0.6 x95 x95 0.6], yl([1 1 2 2]), [0.90 0.90 0.90], 'EdgeColor','none', 'HandleVisibility','off');
patch(ax, [x95 x77 x77 x95], yl([1 1 2 2]), [0.83 0.88 0.95], 'EdgeColor','none', 'HandleVisibility','off');
patch(ax, [x77 6.4 6.4 x77], yl([1 1 2 2]), [0.72 0.80 0.93], 'EdgeColor','none', 'HandleVisibility','off');

for q = 1:numel(machines)
    M = T(T.snr_db == machines(q), :);
    jit = (q-2.5)*0.06;                                  % small offset so the dots do not overlap
    scatter(ax, M.lvl + jit, M.margin_2frame_p95_dB, 8, col(q,:), 'filled', ...
        'MarkerFaceAlpha', 0.25, 'HandleVisibility','off');
    med = arrayfun(@(k) median(M.margin_2frame_p95_dB(M.lvl == k)), 1:numel(levels));
    plot(ax, (1:numel(levels)) + jit, med, '-o', 'Color', col(q,:), 'MarkerFaceColor', col(q,:), ...
        'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', names(q) + ", median");
end
yline(ax, 0, ':', 'Color', [0.80 0.20 0.20], 'LineWidth', 2, 'HandleVisibility','off');   % detection threshold
xticks(ax, 1:numel(levels)); xticklabels(ax, compose('%.1f\\newline(%g%%)', levels/100*pi/6*20^3, levels));
ax.TickLabelInterpreter = 'tex';            % volume (mm^3), share of region underneath
xlim(ax, [0.6 6.4]); ylim(ax, yl);
xlabel(ax, 'Extra fluid in a 20 mm region of brain (mm^3)');
ylabel(ax, 'Margin above the detection threshold (dB)');
legend(ax, 'Location','northwest'); legend(ax,'boxoff');
exportgraphics(f, fullfile(root,'figures','aim2','a2_margin_by_fluid.png'), 'Resolution', 300);
