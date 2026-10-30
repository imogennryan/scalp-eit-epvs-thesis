%A1_VOLUME_LAW_FIGURE  Plot how the change in the readings grows with the radius of a fluid sphere.

% Aim served: Aim 1.
%
% Inputs : results/aim1/a1_analytic_volume_law.csv (from a1_analytic_sphere_validation)
% Outputs: figures/aim1/a1_volume_law.png
% Imogen Ryan - Honours thesis.

clear; close all force;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
T = readtable(fullfile(root,'results','aim1','a1_analytic_volume_law.csv'));

f = figure('Units','centimeters','Position',[1 1 17 13],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
tl = tiledlayout(f, 2, 1, 'Padding','compact', 'TileSpacing','compact');   % stacked so the labels fit at 10 pt

cases = {'low',  0.101, [0.20 0.40 0.70], 'o', '(a) Sphere 0.101 S/m (1% more conductive than brain)';
         'high', 1.79,  [0.85 0.33 0.10], 's', '(b) Sphere 1.79 S/m (CSF, 18 times brain)'};
for c = 1:2
    ax = nexttile(tl); hold(ax,'on')
    set(ax,'XScale','log','YScale','log'); grid(ax,'off'); box(ax,'off')
    m  = startsWith(T.contrast, cases{c,1});
    a  = T.a_mm(m);  s = T.norm_dv(m)*1e6;              % radius (mm), change in readings (uV)
    small = a <= 8;                                      % the law holds for small spheres
    p  = polyfit(log(a(small)), log(s(small)), 1);       % slope on log-log axes
    % physics expectation through the smallest sphere: change grows with radius cubed
    plot(ax, a, s(1)*(a/a(1)).^3, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 2.5, ...
        'DisplayName', 'Physics expectation: slope 3 (change \propto radius^3)');
    plot(ax, a, s, '-', 'Color', cases{c,3}, 'LineWidth', 1.5, 'Marker', cases{c,4}, ...
        'MarkerFaceColor', cases{c,3}, 'MarkerSize', 7, ...
        'DisplayName', sprintf('Model: slope %.3f (fitted 3 to 8 mm)', p(1)));
    xticks(ax, [3 4 5 6 8 10 13 16 20]); xlim(ax, [2.8 21]);
    xlabel(ax, 'Radius of the fluid sphere (mm)');
    ylabel(ax, 'Change in the readings, ||\DeltaV|| (\muV)');
    title(ax, cases{c,5});
    legend(ax, 'Location','northwest'); legend(ax,'boxoff');
end
exportgraphics(f, fullfile(root,'figures','aim1','a1_volume_law.png'), 'Resolution', 300);
