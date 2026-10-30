%A2_NOISE_METHOD_FIGURE  Show how the detection threshold is set from noise alone.

% Aim served: Aim 2 (Methods figure).
%
% Inputs : data/aim2/a2_sweep_m19_1020.mat, results/aim2/a2_sweep_index.csv,
%          a2_noise_model (same seeds as a2_add_noise)
% Outputs: figures/methods/m7_noise_method.png
% Imogen Ryan - Honours thesis.

clear; close all force;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'code','aim2_simulations'));
out  = fullfile(root,'figures','methods'); if ~exist(out,'dir'), mkdir(out); end

% Healthy readings and the fluid case: 10% fluid, 40 mm below Cz, 19 electrodes
S   = load(fullfile(root,'data','aim2','a2_sweep_m19_1020.mat'),'vb');
vb  = S.vb;                                             % healthy readings (V)
T   = readtable(fullfile(root,'results','aim2','a2_sweep_index.csv'),'TextType','string');
key = "B|m19_1020|Cz|d40|s20|b0.1|maxwell_garnett";
pk_fluid = T.peak_dV_mV(T.key == key) * 1e3;
assert(isscalar(pk_fluid), 'expected one matching run');            % fluid's peak change (uV)

% Noise, exactly as in a2_add_noise.m
SNR  = [45 48 77 95];                                   % all four levels, same draw order as a2_add_noise
NA   = a2_noise_model(vb, SNR, 200, 0);                 % noise on the healthy frame
NB   = a2_noise_model(vb, SNR, 200, 1);                 % noise on the second frame
dn   = (NB(:,:,[2 4]) - NA(:,:,[2 4])) * 1e6;           % change from noise alone, 48 and 95 dB (uV)
pk   = squeeze(max(abs(dn), [], 1));                    % largest of the 304, per draw (200 x 2)
bar95 = [prctile_p95(pk(:,1)) prctile_p95(pk(:,2))];    % detection threshold (uV)
sigma = rms(vb)*1e6 / 10^(48/20);                       % noise per reading, 48 dB (uV)

f = figure('Units','centimeters','Position',[1 1 17 12],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
tl = tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
blue = [0.20 0.45 0.75]; red = [0.80 0.20 0.20]; green = [0.20 0.60 0.30];

% (a) One pair of noisy frames, no fluid: the change noise alone gives
ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
d1 = dn(:,1,1); [m1, i1] = max(abs(d1));
stem(ax, 1:numel(d1), d1, 'Marker','none', 'Color',blue);
plot(ax, i1, d1(i1), 'o', 'MarkerFaceColor',red, 'MarkerEdgeColor','k', 'MarkerSize',8);
text(ax, i1+8, d1(i1), sprintf('Largest: %.1f \\muV', m1), 'Color',red, 'FontWeight','bold');
xlim(ax,[0 numel(d1)+1]);
xlabel(ax,'Reading number'); ylabel(ax,'Change from noise alone (\muV)');
title(ax,'(a) Two noise frames, no fluid','FontWeight','normal');

% (b) Repeat 200 times: the threshold is the 95th percentile of the largest values
ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
histogram(ax, pk(:,1), 20, 'FaceColor',blue, 'EdgeColor','w');
xline(ax, bar95(1), '--', 'Color',red, 'LineWidth',2);
text(ax, bar95(1), ax.YLim(2)*0.75, sprintf('Threshold %.1f \\muV ', bar95(1)), ...
     'Color',red, 'FontWeight','bold', 'HorizontalAlignment','right');
xlabel(ax,'Largest change from noise alone (\muV)'); ylabel(ax,'Number of repeats');
title(ax,'(b) Largest change, 200 repeats','FontWeight','normal');

% (c) Compare the fluid's peak change with the threshold: margin in dB
ax = nexttile(tl, 3, [1 2]); hold(ax,'on'); grid(ax,'on'); box(ax,'off');   % full width below (a) and (b)
vals = [bar95(2) pk_fluid bar95(1)];                    % bottom to top: 95 dB threshold, fluid, 48 dB threshold
b = barh(ax, 1:3, vals, 0.55, 'FaceColor','flat', 'EdgeColor','none');
b.CData = [red*0.6+0.4; green; red];
set(ax,'XScale','log','YTick',1:3,'YTickLabel',{'Threshold at 95 dB','Fluid signal','Threshold at 48 dB'});
xlim(ax,[0.03 300]); xticks(ax,[0.1 1 10 100]); xticklabels(ax,{'0.1','1','10','100'}); ylim(ax,[0.4 3.8]);
for k = 1:3, text(ax, vals(k)*1.25, k, sprintf('%.3g \\muV', vals(k)), 'VerticalAlignment','middle'); end
mg = 20*log10(pk_fluid ./ bar95);                       % margin (dB): [48 dB, 95 dB]
% arrows from the fluid signal to each threshold, with the margin written above
plot(ax, [pk_fluid bar95(1)], [2.5 2.5], '-', 'Color',red, 'LineWidth',1.2);
plot(ax, bar95(1), 2.5, '>', 'MarkerFaceColor',red, 'MarkerEdgeColor',red, 'MarkerSize',6);
text(ax, sqrt(pk_fluid*bar95(1)), 2.62, sprintf('%+.0f dB, hidden', mg(1)), 'HorizontalAlignment','center','FontWeight','bold','Color',red);
plot(ax, [bar95(2) pk_fluid], [1.5 1.5], '-', 'Color',green, 'LineWidth',1.2);
plot(ax, bar95(2), 1.5, '<', 'MarkerFaceColor',green, 'MarkerEdgeColor',green, 'MarkerSize',6);
text(ax, sqrt(pk_fluid*bar95(2)), 1.62, sprintf('%+.0f dB, detected', mg(2)), 'HorizontalAlignment','center','FontWeight','bold','Color',green);
xlabel(ax,'Peak change (\muV, log scale)');
title(ax,'(c) Margin = 20 log_{10}(fluid \div threshold)','FontWeight','normal');

exportgraphics(f, fullfile(out,'m7_noise_method.png'),'Resolution',300);
fprintf('Threshold 48 dB %.2f uV, 95 dB %.4f uV. Fluid %.3f uV. Margins %.1f, %.1f dB\n', ...
    bar95, pk_fluid, mg);


function q = prctile_p95(x)
% 95th percentile, same rule as prctile_local in a2_add_noise.m
x = sort(x(:)); n = numel(x);
p = 100*((1:n)' - 0.5)/n;
q = interp1(p, x, 95, 'linear', 'extrap');
end
