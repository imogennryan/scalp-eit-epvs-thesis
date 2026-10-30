%A1_SATURATION_FIGURE  Plot how the change in the readings levels off as a sphere becomes more conductive.

% Aim served: Aim 1.
%
% Inputs : results/aim1/a1_analytic_saturation_law.csv (from a1_analytic_sphere_validation)
% Outputs: figures/aim1/a1_saturation_law.png
% Imogen Ryan - Honours thesis.

clear; close all force;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
T = readtable(fullfile(root,'results','aim1','a1_analytic_saturation_law.csv'));

% Sphere more conductive than brain; the two reference lines pass through the 2x sphere
T = T(T.ratio >= 2, :);
ref = T.ratio == 2;
sb    = T.sigma_i(ref)/2;                          % brain conductivity (S/m)
model = T.resp*1e6;                                % change in the readings (uV)
exact = model(ref) * T.beta/T.beta(ref);           % physics expectation: (si - sb)/(si + 2 sb)
guess = model(ref) * T.naive_linear/T.naive_linear(ref);  % proportional to conductivity: (si - sb)/sb

f = figure('Units','centimeters','Position',[1 1 14 10],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
ax = axes(f); hold(ax,'on')
set(ax,'XScale','log','YScale','log'); grid(ax,'off'); box(ax,'off')
plot(ax, T.sigma_i, guess, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 2, ...
    'DisplayName', 'Signal proportional to conductivity');
plot(ax, T.sigma_i, exact, '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5, ...
    'DisplayName', 'Physics expectation: levels off');
plot(ax, T.sigma_i, model, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, ...
    'Marker', 's', 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 7, ...
    'DisplayName', 'Model');
xticks(ax, [0.2 0.4 0.8 1.6 4 10 40]); xlim(ax, [0.18 45]);
xlabel(ax, 'Sphere conductivity (S/m); brain is 0.1 S/m');
ylabel(ax, 'Change in the readings (\muV)');
legend(ax, 'Location','northwest'); legend(ax,'boxoff');
exportgraphics(f, fullfile(root,'figures','aim1','a1_saturation_law.png'), 'Resolution', 300);
