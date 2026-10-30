%A2_SENSITIVITY_FIGURE  Plot how much the peak change grows from worst to best case for each setting.

% Aim served: Aim 2.
%
% Inputs : results/aim2/a2_results_summary.csv (from a2_results_summary.m),
%          the paired effects with everything else held the same
% Outputs: figures/aim2/a2_sensitivity_worst_best.png
% Imogen Ryan - Honours thesis.

clear; close all force;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
S = readtable(fullfile(root,'results','aim2','a2_results_summary.csv'), 'TextType','string');
get = @(q, sub) S(S.quantity == q & contains(S.subset, sub), :);

rows = {  % label, worst -> best, quantity, subset
 'Extra fluid', '10.5 \rightarrow 418.9 mm^3', 'c: effect of fluid volume', '10% vs 0.25%'
 'Pocket shape', 'round pockets (main case) \rightarrow tubes lined up with the current', 'c: effect of inclusion shape', 'primary geometry (scenario B) | L=0 vs L=1/3'
 'Scalp point', 'T7 \rightarrow Fz (25-35 mm deep)', 'c: effect of position', 'spread across 6 anchors, 25-35 mm'
 'Electrode layout', '16 or 32 electrodes \rightarrow 19 electrodes', 'c: effect of montage (on margin)', '| 48 dB'
 'Depth', 'deepest that fits \rightarrow 25 mm', 'c: effect of depth', 'shallowest (25 mm) vs deepest reachable'
 'Region size, same fluid', 'weakest \rightarrow strongest of five sizes (12 to 32 mm)', 'c/d: effect of region size at fixed fluid volume', 'L=1/3'};
n = size(rows,1); med = zeros(n,1); sd = med;
for k = 1:n
    r = get(rows{k,3}, rows{k,4}); med(k) = r.median(1); sd(k) = r.sd(1);
end

f = figure('Units','centimeters','Position',[1 1 17 10],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
ax = axes(f); hold(ax,'on'); grid(ax,'off'); box(ax,'off')
y = (n:-1:1) * 1.4;                                  % room above each bar for its label
barh(ax, y, med, 0.45, 'FaceColor', [0.30 0.50 0.75], 'EdgeColor','none');
errorbar(ax, med, y, sd, 'horizontal', 'k', 'LineStyle','none', 'CapSize', 6);
for k = 1:n
    text(ax, 0.3, y(k) + 0.5, ['\bf' rows{k,1} ':\rm  ' rows{k,2}], 'FontSize', 11, 'Interpreter','tex');
    text(ax, med(k) + sd(k) + 0.5, y(k), sprintf('%.1f dB  (%.3g\\times)', med(k), 10^(med(k)/20)), ...
        'VerticalAlignment','middle', 'FontSize', 11, 'FontWeight','bold');
end
yticks(ax, []); ylim(ax, [0.6, y(1) + 0.9]); xlim(ax, [0 40]);
xlabel(ax, 'Increase in peak change from worst case to best case (dB); median with 1 SD');
exportgraphics(f, fullfile(root,'figures','aim2','a2_sensitivity_worst_best.png'), 'Resolution', 300);
disp(table(string(rows(:,1)), med, sd, 10.^(med/20), 'VariableNames', {'setting','median_dB','sd_dB','ratio'}))
