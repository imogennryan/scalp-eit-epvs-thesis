%A2_RESULTS_FIGURES_V3  Draw the Aim 2 Results figures from saved outputs.

% Aim served: Aim 2.
%
% Reads saved results only, with no forward solves. Main case: a 20 mm region
% holding 0.25 to 10% fluid, Maxwell Garnett mixing (L = 1/3). The margin
% compares the peak change with the 95th percentile of the peak from noise
% alone over two frames. The map shows detection only, not location.
% Noise levels: clinical 45 and 48 dB (Goren 2018, McDermott 2020), research
% 77 and 95 dB (Ouypornkochagorn 2022).
%
% Inputs : results/aim2/a2_sweep_index.csv, a2_noise_summary.csv, a2_noise_floor.csv,
%          a2_sweep_reachability.csv, a2_results_summary.csv,
%          data/aim2/a2_sweep_m19_1020.mat
% Outputs: figures/aim2/a2_detectability_map_v3 (.png and .fig),
%          results/aim2/a2_figure_numbers_v3.csv
% Imogen Ryan - Honours thesis.

clear; close all force;
rng(0, 'twister');   % fixed seed, although no random draws are made

root = thesis_root();
resd = fullfile(root, 'results', 'aim2');
datd = fullfile(root, 'data', 'aim2');
figd = fullfile(root, 'figures', 'aim2');

PRIMARY  = "maxwell_garnett";                    % L = 1/3
SNR_DB   = [45 48 77 95];                        % instrument SNR (dB), see header
SNR_NAME = ["clinical, higher noise", "clinical, lower noise", "research, higher noise", "research, lowest noise"];
BURDENS  = [0.0025 0.005 0.01 0.02 0.05 0.10];   % fluid share of the region, as in a2_run_sweep
MONTAGES = ["m19_1020", "m32_1010", "m16_twoplane"];

F1 = struct('montage', "m19_1020", 'anchor', "Cz", 'depth_mm', 40, ...
            'burdens', [0.0025 0.01 0.02 0.05 0.10], 'example_burden', 0.02);
% Research instrument on top
F3 = struct('montage', "m19_1020", 'burden', 0.10, 'snr', [77 48], ...
            'anchors', ["Fz","Cz","C3","Pz","T7","O1"], 'depths', 25:5:55);

C.black = [0 0 0];

fprintf('\n########### A2: RESULTS FIGURES v3 (reads saved outputs only) ###########\n');

NUM = table('Size', [0 5], 'VariableTypes', ...
    {'string','string','string','double','string'}, 'VariableNames', ...
    {'figure','quantity','subset','value','units'});
addnum = @(T, f, q, s, v, u) [T; table(string(f), string(q), string(s), ...
    double(v), string(u), 'VariableNames', T.Properties.VariableNames)];

%% Load (read only)
I   = readtable(fullfile(resd, 'a2_sweep_index.csv'),        'TextType', 'string');
N   = readtable(fullfile(resd, 'a2_noise_summary.csv'),      'TextType', 'string');
F   = readtable(fullfile(resd, 'a2_noise_floor.csv'),        'TextType', 'string');
RCH = readtable(fullfile(resd, 'a2_sweep_reachability.csv'), 'TextType', 'string');
RS  = readtable(fullfile(resd, 'a2_results_summary.csv'),    'TextType', 'string');
S = load(fullfile(datd, char("a2_sweep_" + F1.montage + ".mat")), 'DV', 'vb');
vb = S.vb; DV = S.DV;

PB = I(I.scenario == "B" & I.rule == PRIMARY, :);
P  = N(N.scenario == "B" & N.rule == PRIMARY, :);

%% Integrity checks, before any figure
fprintf('\n=============== INTEGRITY CHECKS ===============\n');
ok = true;
ok = chk(ok, 'index rows', height(I), height(I) == 1578, 'expected 1578');
ok = chk(ok, 'primary (scenario B, L=1/3) rows', height(PB), height(PB) == 558, ...
    'expected 558 = 186 per montage x 3');
ok = chk(ok, 'noise summary rows', height(N), height(N) == 4*height(I), ...
    'expected 4 SNR levels per configuration');
nb = sum(~isfinite(PB.peak_dV_mV)) + sum(PB.peak_dV_mV == 0) + ...
     sum(~isfinite(P.margin_2frame_p95_dB)) + sum(~isfinite(PB.fluid_vol_mm3));
ok = chk(ok, 'NaN / Inf / exact-zero in peak, margin, fluid', nb, nb == 0, 'primary subset');
ok = chk(ok, 'unique primary keys', numel(unique(PB.key)), ...
    numel(unique(PB.key)) == height(PB), 'duplicates would show as a smaller count');
p95 = nan(height(P), 1);
for j = 1:height(F)
    k = P.montage == F.montage(j) & P.snr_db == F.snr_db(j);
    p95(k) = F.peak_noise_p95_2frame_mV(j);
end
e = max(abs(20*log10(P.peak_dV_mV ./ p95) - P.margin_2frame_p95_dB));
ok = chk(ok, 'margin re-derived from a2_noise_floor.csv (dB)', e, e < 1e-9, ...
    'margin = 20 log10(peak / p95 two-frame noise-only peak)');
sel = find(PB.montage == F1.montage);
pk = zeros(numel(sel), 1);
for q = 1:numel(sel), pk(q) = max(abs(DV(char(PB.key(sel(q)))))) * 1e3; end
e = max(abs(pk - PB.peak_dV_mV(sel)) ./ PB.peak_dV_mV(sel));
ok = chk(ok, 'peak from stored vectors vs index (rel.)', e, e < 1e-9, char(F1.montage));
nmis = 0;
for a = F3.anchors
    dmax = RCH.max_centre_region_mm(RCH.anchor == a);
    for d = F3.depths
        have = any(P.montage == F3.montage & P.anchor == a & P.depth_mm == d & ...
                   abs(P.burden_frac - F3.burden) < 1e-12 & P.snr_db == F3.snr(1));
        nmis = nmis + (have ~= (d <= dmax));
    end
end
ok = chk(ok, 'map cells agree with a2_sweep_reachability.csv', nmis, nmis == 0, ...
    'a cell exists iff a 20 mm region centred there fits inside the brain');
nclipped = sum(PB.fits_in_brain == 0);
fprintf('  NOTE     %d of %d primary rows delivered < 99%% of the nominal region volume\n', ...
    nclipped, height(PB));
fprintf('           (T7 and O1 at 25 mm, 95.5-97.7%%; kept, within the 5%% tolerance of a2_check_dataset.m)\n');
fprintf('\nALL CHECKS PASS: %d\n', ok);
assert(ok, 'An integrity check failed. Do not use any figure below until it is understood.');

%% Figure 4.1: healthy against diseased, change per reading
fprintf('\n=============== FIGURE 4.1: HEALTHY vs DISEASED (uV) ===============\n');
nb = numel(F1.burdens);
M = zeros(nb, 3); SD = zeros(nb, 3); fluid = zeros(nb, 1);
for q = 1:nb
    k = PB.montage == F1.montage & PB.anchor == F1.anchor & ...
        PB.depth_mm == F1.depth_mm & abs(PB.burden_frac - F1.burdens(q)) < 1e-12;
    dv = DV(char(PB.key(k)));                % difference vector (V)
    v1 = vb + dv;                            % readings with fluid (V)
    fluid(q) = PB.fluid_vol_mm3(k);
    M(q, :)  = [mean(abs(vb)) mean(abs(v1)) mean(abs(dv))] * 1e6;   % uV
    SD(q, :) = [std(abs(vb))  std(abs(v1))  std(abs(dv))]  * 1e6;   % uV
    fprintf('  %5.2f %% v/v | %6.1f mm3 | healthy %.2f+/-%.2f uV | fluid %.2f+/-%.2f uV | change %.4f+/-%.4f uV\n', ...
        F1.burdens(q)*100, fluid(q), M(q,1), SD(q,1), M(q,2), SD(q,2), M(q,3), SD(q,3));
    sub = sprintf('%s | %s | %d mm | %g %% v/v', F1.montage, F1.anchor, F1.depth_mm, F1.burdens(q)*100);
    NUM = addnum(NUM, 'Fig 4.1a', 'fluid volume in the region', sub, fluid(q), 'mm^3');
    lbl = ["mean |healthy reading|", "mean |reading with fluid|", "mean |change|"];
    for c = 1:3
        NUM = addnum(NUM, 'Fig 4.1a', lbl(c), sub, M(q, c), 'uV');
        NUM = addnum(NUM, 'Fig 4.1a', lbl(c) + ", SD over readings", sub, SD(q, c), 'uV');
    end
end

% Panel (b): per-reading change against the noise-alone peaks
ke = PB.montage == F1.montage & PB.anchor == F1.anchor & ...
     PB.depth_mm == F1.depth_mm & abs(PB.burden_frac - F1.example_burden) < 1e-12;
ex = PB(ke, :);
adv = abs(DV(char(ex.key))) * 1e6;                        % uV
fl = @(s) 1e3 * F.peak_noise_p95_2frame_mV(F.montage == F1.montage & F.snr_db == s);
f48 = fl(48); f95 = fl(95);
above = adv > f95;
[pkv, ipk] = max(adv);
assert(nnz(adv > f48) == 0, 'A reading now clears the clinical noise level; revisit panel b.');

fprintf('  (b) %s: peak %.4f uV; noise-alone p95 48 dB %.4f uV, 95 dB %.4f uV; %d of %d above 95 dB\n', ...
    ex.key, pkv, f48, f95, nnz(above), numel(adv));
NUM = addnum(NUM, 'Fig 4.1b', 'peak change', char(ex.key), pkv, 'uV');
NUM = addnum(NUM, 'Fig 4.1b', 'noise-alone peak, p95, 48 dB', char(F1.montage), f48, 'uV');
NUM = addnum(NUM, 'Fig 4.1b', 'noise-alone peak, p95, 95 dB', char(F1.montage), f95, 'uV');
NUM = addnum(NUM, 'Fig 4.1b', 'readings above the 95 dB level', char(ex.key), nnz(above), 'count');
NUM = addnum(NUM, 'Fig 4.1b', 'readings in total', char(ex.key), numel(adv), 'count');

% Panel (c): least-squares line of log10 peak on log10 fluid volume, as in a2_results_summary
x  = log10(PB.fluid_vol_mm3); y = log10(PB.peak_dV_mV * 1e3);   % peak in uV
n  = numel(x);
bc = [ones(n, 1) x] \ y;
res = y - [ones(n, 1) x]*bc;
r2 = 1 - sum(res.^2) / sum((y - mean(y)).^2);
se = sqrt(sum(res.^2) / (n - 2) / sum((x - mean(x)).^2));
kS = RS.quantity == "e: log-log slope, peak |dV| vs fluid volume" & RS.subset == "primary | pooled";
kR = RS.quantity == "e: R^2, variance explained by fluid volume"   & RS.subset == "primary | pooled";
fprintf('  (c) pooled OLS: slope %.4f (SE %.4f), R^2 %.4f, r %.4f, n %d\n', bc(2), se, r2, sqrt(r2), n);
fprintf('      a2_results_summary.csv: slope %.4f, R^2 %.4f, n %d\n', RS.median(kS), RS.median(kR), RS.n(kS));
assert(abs(bc(2) - RS.median(kS)) < 1e-9 && abs(r2 - RS.median(kR)) < 1e-9 && n == RS.n(kS), ...
    'Recomputed slope / R^2 / n do not match a2_results_summary.csv. Stop and investigate.');
for m = MONTAGES
    k = PB.montage == m;
    xm = log10(PB.fluid_vol_mm3(k)); ym = log10(PB.peak_dV_mV(k) * 1e3);
    bm = [ones(nnz(k), 1) xm] \ ym; rm = ym - [ones(nnz(k), 1) xm]*bm;
    r2m = 1 - sum(rm.^2) / sum((ym - mean(ym)).^2);
    fprintf('      %-13s slope %.4f  R^2 %.4f  r %.4f  n %d\n', m, bm(2), r2m, sqrt(r2m), nnz(k));
    NUM = addnum(NUM, 'Fig 4.1c', 'log-log slope', char(m), bm(2), '-');
    NUM = addnum(NUM, 'Fig 4.1c', 'R^2', char(m), r2m, 'fraction');
end
NUM = addnum(NUM, 'Fig 4.1c', 'log-log slope', 'pooled', bc(2), '-');
NUM = addnum(NUM, 'Fig 4.1c', 'log-log slope, SE', 'pooled', se, '-');
NUM = addnum(NUM, 'Fig 4.1c', 'R^2', 'pooled', r2, 'fraction');
NUM = addnum(NUM, 'Fig 4.1c', 'configurations', 'pooled', n, 'count');

%% Figure 4.2: margin against fluid volume
fprintf('\n=============== FIGURE 4.2: MARGIN vs FLUID VOLUME ===============\n');
XB = zeros(numel(BURDENS), numel(SNR_DB)); YB = XB;
for is = 1:numel(SNR_DB)
    k = P.snr_db == SNR_DB(is);
    [~, loc] = ismember(P.key(k), PB.key);
    fv = PB.fluid_vol_mm3(loc); mg = P.margin_2frame_p95_dB(k);
    for q = 1:numel(BURDENS)
        kb = abs(P.burden_frac(k) - BURDENS(q)) < 1e-12;
        XB(q, is) = median(fv(kb)); YB(q, is) = median(mg(kb));
    end
end
xc = nan(1, numel(SNR_DB));
for is = 1:numel(SNR_DB)
    j = find(YB(1:end-1, is) < 0 & YB(2:end, is) >= 0, 1);
    if ~isempty(j)
        t = (0 - YB(j, is)) / (YB(j+1, is) - YB(j, is));
        xc(is) = 10^(log10(XB(j, is)) + t*(log10(XB(j+1, is)) - log10(XB(j, is))));
    end
    fprintf('  %2d dB: median crosses 0 dB at %.1f mm3 (NaN = never)\n', SNR_DB(is), xc(is));
    for q = 1:numel(BURDENS)
        NUM = addnum(NUM, 'Fig 4.2', 'median margin above noise alone', ...
            sprintf('%d dB | fluid %.1f mm^3', SNR_DB(is), XB(q, is)), YB(q, is), 'dB');
    end
    NUM = addnum(NUM, 'Fig 4.2', 'fluid volume where the median crosses 0 dB', ...
        sprintf('%d dB (%s)', SNR_DB(is), SNR_NAME(is)), xc(is), 'mm^3');
end
assert(all(isnan(xc([1 2]))) && all(~isnan(xc([3 4]))), 'Band layout assumption broken.');

%% Figure 4.3: detectability map, 77 dB on top
fprintf('\n=============== FIGURE 4.3: DETECTABILITY MAP ===============\n');
na = numel(F3.anchors); nd = numel(F3.depths);
MAP = nan(na, nd, numel(F3.snr));
for is = 1:numel(F3.snr)
    for ia = 1:na
        for id = 1:nd
            k = P.montage == F3.montage & P.anchor == F3.anchors(ia) & ...
                P.depth_mm == F3.depths(id) & abs(P.burden_frac - F3.burden) < 1e-12 & ...
                P.snr_db == F3.snr(is);
            if any(k), MAP(ia, id, is) = P.margin_2frame_p95_dB(k); end
        end
    end
end
clim_max = ceil(max(abs(MAP(:)), [], 'omitnan') / 5) * 5;
cmap = red_green_map(256);

f3 = newfig([4 2 17 17]);                                % printed at 17 cm, so 10 pt here is 10 pt on the page
tl = tiledlayout(f3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
pan = 'ab';
for is = 1:numel(F3.snr)
    ax = nexttile(tl); hold(ax, 'on');
    Z = MAP(:, :, is);
    imagesc(ax, 1:nd, 1:na, Z, 'AlphaData', ~isnan(Z));
    colormap(ax, cmap); clim(ax, [-clim_max clim_max]);
    for ia = 1:na
        for id = 1:nd
            if isnan(Z(ia, id))
                text(ax, id, ia, sprintf('Outside\nbrain'), 'FontSize', 10, ...
                    'Color', [0.55 0.55 0.55], 'HorizontalAlignment', 'center');
            else
                tc = C.black;
                if abs(Z(ia, id)) > 0.55*clim_max, tc = [1 1 1]; end
                text(ax, id, ia, sprintf('%+.1f', Z(ia, id)), 'FontSize', 10, ...
                    'Color', tc, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            end
        end
    end
    for id = 0.5:1:(nd + 0.5)
        plot(ax, [id id], [0.5 na + 0.5], '-', 'Color', [1 1 1], 'LineWidth', 1.5);
    end
    for ia = 0.5:1:(na + 0.5)
        plot(ax, [0.5 nd + 0.5], [ia ia], '-', 'Color', [1 1 1], 'LineWidth', 1.5);
    end
    set(ax, 'YDir', 'reverse', 'XTick', 1:nd, 'XTickLabel', compose('%d', F3.depths'), ...
        'YTick', 1:na, 'YTickLabel', cellstr(F3.anchors'), 'TickLength', [0 0]);
    xlim(ax, [0.5 nd + 0.5]); ylim(ax, [0.5 na + 0.5]);
    ylabel(ax, 'Scalp point');
    if is == numel(F3.snr), xlabel(ax, 'Depth of the region below the scalp point (mm)'); end
    ndet = nnz(Z > 0); ntot = nnz(~isnan(Z));
    mname = containers.Map({45, 48, 77, 95}, {'Clinical, higher noise', 'Clinical, lower noise', ...
        'Research, higher noise', 'Research, lowest noise'});             % machine names used in the report
    title(ax, sprintf('(%s)  %s (%d dB)', pan(is), mname(F3.snr(is)), F3.snr(is)), 'FontWeight', 'bold');
    ax.TitleHorizontalAlignment = 'left';
    figstyle(ax, 10);
    set(ax, 'Box', 'on');
    fprintf('  (%s) %d dB: %d of %d reachable positions above 0 dB; range %+.1f to %+.1f dB\n', ...
        pan(is), F3.snr(is), ndet, ntot, min(Z(:), [], 'omitnan'), max(Z(:), [], 'omitnan'));
    sub = sprintf('%s | %g %% v/v | %d dB', F3.montage, F3.burden*100, F3.snr(is));
    NUM = addnum(NUM, 'Fig 4.3', 'reachable positions above 0 dB', sub, ndet, 'count');
    NUM = addnum(NUM, 'Fig 4.3', 'reachable positions in total', sub, ntot, 'count');
    NUM = addnum(NUM, 'Fig 4.3', 'weakest margin', sub, min(Z(:), [], 'omitnan'), 'dB');
    NUM = addnum(NUM, 'Fig 4.3', 'strongest margin', sub, max(Z(:), [], 'omitnan'), 'dB');
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
cb.Label.String = 'Margin above noise alone (dB)';
cb.FontSize = 10; cb.Label.FontSize = 10;
cb.Ticks = -clim_max:10:clim_max;
NUM = addnum(NUM, 'Fig 4.3', 'colour scale limit (symmetric)', 'both panels', clim_max, 'dB');
savefig_both(f3, fullfile(figd, 'a2_detectability_map_v3'));

%% Numbers behind the figures, and the re-run fingerprint
numcsv = fullfile(resd, 'a2_figure_numbers_v3.csv');
writetable(NUM, numcsv);
Nb = readtable(numcsv, 'TextType', 'string');
assert(height(Nb) == height(NUM) && max(abs(Nb.value - NUM.value), [], 'omitnan') < 1e-9, ...
    'Figure-numbers CSV read-back does not match.');
fprintf('\n  %d numbers written to %s\n  FINGERPRINT: %010d\n', height(NUM), numcsv, fingerprint(NUM.value));
fprintf('\nSAVED figures/aim2/a2_detectability_map_v3.{png,fig}\n');

%% Local functions
function ok = chk(ok, name, value, pass, detail)
if pass, v = 'PASS'; else, v = '**FAIL**'; end
fprintf('  %-8s %-50s value %-12.4g %s\n', v, name, double(value), detail);
ok = ok && pass;
end

function h = fingerprint(v)
s = sprintf('%.12g|', v(:));
h = uint64(5381);
for j = 1:numel(s)
    h = mod(h*uint64(33) + uint64(double(s(j))), uint64(2147483647));
end
h = double(h);
end

function s = capitalise(x)
s = char(x); s(1) = upper(s(1));
end

function cmap = red_green_map(n)
% Red below noise, pale at 0 dB, green above noise
lo  = [0.78 0.12 0.12];      % red
mid = [0.96 0.95 0.90];      % pale neutral
hi  = [0.10 0.55 0.22];      % green
m = floor(n/2); t = linspace(0, 1, m)';
cmap = [lo + t.*(mid - lo); mid + t.*(hi - mid)];
end

function f = newfig(pos_cm)
f = figure('Color', 'w', 'Units', 'centimeters', 'Position', pos_cm);
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
theme(f, 'light');
end

function figstyle(ax, fs)
set(ax, 'FontSize', fs, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1, ...
    'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], 'Layer', 'top');
grid(ax, 'off');
end

function savefig_both(fig, base)
exportgraphics(fig, [base '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
savefig(fig, [base '.fig']);
end

function r = thesis_root()
r = fileparts(mfilename('fullpath'));
while ~exist(fullfile(r, 'README.md'), 'file')
    p = fileparts(r);
    if isempty(p) || strcmp(p, r), error('thesis_root:noRoot', 'No README.md found.'); end
    r = p;
end
end
