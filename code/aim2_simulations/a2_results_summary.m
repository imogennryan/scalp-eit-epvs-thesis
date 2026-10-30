%A2_RESULTS_SUMMARY  Summarise the simulated results, healthy head against head with fluid, with an SD beside every median.

% Aim served: Aim 2.
%
% Reads saved outputs only; no forward solve. Integrity checks run first.
% Control is the healthy head (V0); experimental is the same head with one
% region of extra fluid (V1 = V0 + dV), on the identical mesh and electrodes.
% Main case: scenario B, maxwell_garnett (L = 1/3), two-frame noise, margin
% against the 95th percentile of the noise-only peak, all three montages.
% Scenario A is a bound and is never pooled into a main-case number.
% Each axis effect compares configurations that differ in one factor only.
% Each row also gives skewness G1 and the bimodality coefficient BC, flagged
% when |G1| > 1 or BC > 5/9, where the IQR describes the spread better than SD.
%
% Inputs : results/aim2/a2_sweep_index.csv, a2_noise_summary.csv,
%          a2_noise_floor.csv, a2_size_sweep.csv,
%          data/aim2/a2_sweep_<montage>.mat,
%          models/baseline/a1_baseline_<montage>.mat (prov struct only)
% Outputs: results/aim2/a2_results_summary.csv  one row per reported quantity
% Imogen Ryan - Honours thesis.

clear; close all force;
rng(0, 'twister');   % fixed seed; no random numbers are drawn

root = thesis_root();
resd = fullfile(root, 'results', 'aim2');
datd = fullfile(root, 'data', 'aim2');
basd = fullfile(root, 'models', 'baseline');
figd = fullfile(root, 'figures', 'aim2');
outcsv = fullfile(resd, 'a2_results_summary.csv');

MONTAGES  = ["m19_1020", "m32_1010", "m16_twoplane"];
SNR_DB    = [45 48 77 95];                     % instrument noise levels (dB), Goren 2018, McDermott 2020, Ouypornkochagorn 2022
SNR_NAME  = ["clinical", "clinical", "research", "best research"];
BURDENS   = [0.0025 0.005 0.01 0.02 0.05 0.10];% fluid volume fractions, as in a2_run_sweep
PRIMARY   = "maxwell_garnett";                 % L = 1/3
ALTRULE   = "linear";                          % L = 0
CORE_MAX_DEPTH = 35;                           % deepest depth all 6 anchors reach (mm)
N_EXPECT  = struct('index', 1578, 'B', 1116, 'A', 462, 'per_montage', 526, ...
                   'primary_per_montage', 186, 'size', 180);
% Example configuration for the healthy versus fluid figure. Chosen as typical,
% not best case: Cz, 40 mm deep, 2 % fluid. Its margin percentile is reported below.
EX = struct('montage', "m19_1020", 'anchor', "Cz", 'depth_mm', 40, 'burden_frac', 0.02);

fprintf('\n############### A2: RESULTS SUMMARY (reads saved outputs only) ###############\n');

R = table();   % the output table, built row by row

%% 1. Load
I = readtable(fullfile(resd, 'a2_sweep_index.csv'),   'TextType', 'string');
N = readtable(fullfile(resd, 'a2_noise_summary.csv'), 'TextType', 'string');
F = readtable(fullfile(resd, 'a2_noise_floor.csv'),   'TextType', 'string');
Z = readtable(fullfile(resd, 'a2_size_sweep.csv'),    'TextType', 'string');

VB = struct(); DVM = struct(); PROV = struct();
for m = MONTAGES
    S = load(fullfile(datd, char("a2_sweep_" + m + ".mat")), 'DV', 'vb');
    B = load(fullfile(basd, char("a1_baseline_" + m + ".mat")), 'prov');
    VB.(m) = S.vb;  DVM.(m) = S.DV;  PROV.(m) = B.prov;
end

%% 2. Integrity checks
fprintf('\n=============== INTEGRITY CHECKS ===============\n');
allok = true;

% 2a. The anchor_note column, which contains commas, read as one column.
ok = width(I) == 28 && all(ismember(unique(I.anchor), ["Cz","Fz","C3","Pz","T7","O1"])) ...
     && isnumeric(I.depth_actual_mm) && isnumeric(I.peak_dV_mV);
[R, allok] = addcheck(R, allok, 'index parsed (28 columns, anchor_note intact)', width(I), ok, ...
    sprintf('expected 28 columns, got %d', width(I)));

% 2b. Expected versus actual run counts.
[R, allok] = addcheck(R, allok, 'index rows', height(I), height(I) == N_EXPECT.index, ...
    sprintf('expected %d', N_EXPECT.index));
[R, allok] = addcheck(R, allok, 'scenario B rows', sum(I.scenario == "B"), ...
    sum(I.scenario == "B") == N_EXPECT.B, sprintf('expected %d', N_EXPECT.B));
[R, allok] = addcheck(R, allok, 'scenario A rows', sum(I.scenario == "A"), ...
    sum(I.scenario == "A") == N_EXPECT.A, sprintf('expected %d', N_EXPECT.A));
for m = MONTAGES
    nm = sum(I.montage == m);
    [R, allok] = addcheck(R, allok, char("rows, " + m), nm, nm == N_EXPECT.per_montage, ...
        sprintf('expected %d', N_EXPECT.per_montage));
    np = sum(I.montage == m & I.scenario == "B" & I.rule == PRIMARY);
    [R, allok] = addcheck(R, allok, char("primary rows, " + m), np, ...
        np == N_EXPECT.primary_per_montage, sprintf('expected %d', N_EXPECT.primary_per_montage));
end
[R, allok] = addcheck(R, allok, 'unique index keys', numel(unique(I.key)), ...
    numel(unique(I.key)) == height(I), 'duplicates would show as fewer than the row count');

% 2c. NaN, Inf and exact zeros in the quantities analysed.
bad = @(x) sum(~isfinite(x)) + sum(x == 0);
nb = bad(I.peak_dV_mV) + bad(N.peak_dV_mV) + sum(~isfinite(N.margin_2frame_p95_dB)) + bad(Z.peak_dV_mV);
[R, allok] = addcheck(R, allok, 'NaN / Inf / exact-zero values', nb, nb == 0, ...
    'peak_dV (index, noise summary, size sweep) and margin_2frame_p95_dB');

% 2d. Noise summary is complete and matches the index.
[R, allok] = addcheck(R, allok, 'noise summary rows', height(N), height(N) == 4*height(I), ...
    sprintf('expected %d = %d configurations x 4 SNR', 4*height(I), height(I)));
ukN = unique(N.key + "|" + string(N.snr_db));
ok = numel(ukN) == height(N) && isempty(setxor(unique(N.key), I.key));
[R, allok] = addcheck(R, allok, 'noise summary keys match index, no duplicates', numel(ukN), ok, ...
    'every (key, SNR) once; key set identical to the index');
[tf, loc] = ismember(N.key, I.key);
e = max(abs(N.peak_dV_mV(tf) - I.peak_dV_mV(loc(tf))) ./ I.peak_dV_mV(loc(tf)));
[R, allok] = addcheck(R, allok, 'noise-summary peak vs index peak (max rel. diff)', e, e < 1e-9, ...
    'peak recomputed from vectors in a2_add_noise vs stored in the sweep index');

% 2e. The margin column can be re-derived from the saved noise comparators.
p95 = nan(height(N), 1);
for j = 1:height(F)
    k = N.montage == F.montage(j) & N.snr_db == F.snr_db(j);
    p95(k) = F.peak_noise_p95_2frame_mV(j);
end
e = max(abs(20*log10(N.peak_dV_mV ./ p95) - N.margin_2frame_p95_dB));
[R, allok] = addcheck(R, allok, 'margin re-derived from a2_noise_floor.csv (max abs diff, dB)', e, ...
    e < 1e-9, 'margin = 20 log10(peak / p95 two-frame noise-only peak)');

% 2f. Stored vectors, baselines and provenance.
sig0 = PROV.(MONTAGES(1)).sigma;
for m = MONTAGES
    DV = DVM.(m); vb = VB.(m); sel = find(I.montage == m);
    ok = DV.Count == numel(sel) && isempty(setxor(string(DV.keys)', I.key(sel)));
    [R, allok] = addcheck(R, allok, char("DV map keys = index keys, " + m), DV.Count, ok, ...
        sprintf('expected %d', numel(sel)));
    pk = zeros(numel(sel), 1); nf = 0;
    for q = 1:numel(sel)
        dv = DV(char(I.key(sel(q))));
        nf = nf + sum(~isfinite(dv));
        pk(q) = max(abs(dv)) * 1e3;
    end
    e = max(abs(pk - I.peak_dV_mV(sel)) ./ I.peak_dV_mV(sel));
    [R, allok] = addcheck(R, allok, char("peak from stored vectors vs index, " + m), e, ...
        e < 1e-9 && nf == 0, sprintf('max rel. diff; %d non-finite vector entries', nf));
    ok = numel(vb) == PROV.(m).checks.n_meas && all(I.n_meas(sel) == numel(vb));
    [R, allok] = addcheck(R, allok, char("healthy V0 length = prov n_meas, " + m), numel(vb), ok, ...
        sprintf('prov.checks.n_meas = %d', PROV.(m).checks.n_meas));
    e = max(abs(rms(vb)*1e3 - I.baseline_rms_mV(sel)));
    [R, allok] = addcheck(R, allok, char("RMS(V0) vs index baseline_rms_mV, " + m), e, e < 1e-9, ...
        'max abs diff, mV');
    s = PROV.(m).sigma;
    ok = isequal(s, sig0);
    [R, allok] = addcheck(R, allok, char("baseline conductivity set, " + m), 4, ok, ...
        sprintf('scalp %.2f skull %.2f CSF %.2f brain %.2f S/m; prov generated %s', ...
        s.scalp, s.skull, s.csf, s.brain, PROV.(m).generated));
end

% 2g. Mixing-rule twins: every L = 1/3 run has an L = 0 twin.
twin = replace(I.key(I.rule == PRIMARY), PRIMARY, ALTRULE);
ok = all(ismember(twin, I.key(I.rule == ALTRULE)));
[R, allok] = addcheck(R, allok, 'every L = 1/3 run has an L = 0 twin', numel(twin), ok, 'key match');

% 2h. Size sweep.
[G, ~, ~, ~, ~] = findgroups(Z.anchor, Z.depth_target_mm, Z.fluid_vol_target_mm3, round(Z.L*6));
nd = splitapply(@numel, Z.diam_mm, G);
e = max(abs(Z.fluid_vol_actual_mm3 - Z.fluid_vol_target_mm3) ./ Z.fluid_vol_target_mm3);
ok = height(Z) == N_EXPECT.size && all(nd == 5) && e < 1e-6;
[R, allok] = addcheck(R, allok, 'size sweep complete and volume-matched', height(Z), ok, ...
    sprintf('expected %d rows, %d groups x 5 diameters; worst volume mismatch %.1e', ...
    N_EXPECT.size, max(G), e));

fprintf('\nALL CHECKS PASS: %d\n', allok);
assert(allok, 'An integrity check failed. Do not use any number below until it is understood.');

%% 3. Main-case subsets
PB  = I(I.scenario == "B" & I.rule == PRIMARY, :);   % noiseless, primary
PBL = I(I.scenario == "B" & I.rule == ALTRULE, :);
P   = N(N.scenario == "B" & N.rule == PRIMARY, :);   % with margins, primary
bfmt = @(b) sprintf('%g%%', b*100);

%% (a) Margin above noise alone (dB)
fprintf('\n=============== (a) MARGIN ABOVE NOISE ALONE (dB) ===============\n');
fprintf('Primary scenario. margin = 20 log10(peak|dV| / p95 noise-only peak), two-frame.\n');
for is = 1:numel(SNR_DB)
    s = SNR_DB(is);
    k = P.snr_db == s;
    R = addrow(R, 'a: margin above noise alone', ...
        sprintf('primary | pooled 3 montages | %d dB (%s)', s, SNR_NAME(is)), ...
        P.margin_2frame_p95_dB(k), 'dB', ...
        ['max = best case: ' bestcase(P(k,:))]);
    for m = MONTAGES
        km = k & P.montage == m;
        R = addrow(R, 'a: margin above noise alone', ...
            sprintf('primary | %s | %d dB (%s)', m, s, SNR_NAME(is)), ...
            P.margin_2frame_p95_dB(km), 'dB', ['max = best case: ' bestcase(P(km,:))]);
    end
end
fprintf('\nConditional slice: margin by burden (pooled over montages, anchors, depths)\n');
for is = 1:numel(SNR_DB)
    for b = BURDENS
        k = P.snr_db == SNR_DB(is) & abs(P.burden_frac - b) < 1e-12;
        R = addrow(R, 'a: margin above noise alone, by burden', ...
            sprintf('primary | pooled 3 montages | %d dB | burden %s', SNR_DB(is), bfmt(b)), ...
            P.margin_2frame_p95_dB(k), 'dB', ...
            'conditional slice: spread is anchor + depth + montage at one burden');
    end
end

%% (b) Peak |dV| against the healthy baseline
fprintf('\n=============== (b) PEAK |dV| FROM THE HEALTHY BASELINE (uV) ===============\n');
for m = MONTAGES
    vb = VB.(m);
    R = addrow(R, 'b: healthy baseline |V0| per reading', char("control head | " + m), ...
        abs(vb)*1e6, 'uV', sprintf('spread across the %d healthy readings; RMS = %.1f uV', ...
        numel(vb), rms(vb)*1e6));
    R = addval(R, 'b: healthy baseline RMS', char("control head | " + m), numel(vb), ...
        rms(vb)*1e6, 'uV', 'single value (RMS over all readings) stored in the median column');
end
for b = BURDENS
    k = abs(PB.burden_frac - b) < 1e-12;
    R = addrow(R, 'b: peak |dV| (fluid minus healthy)', ...
        sprintf('primary | pooled 3 montages | burden %s | fluid %.1f mm3', bfmt(b), ...
        median(PB.fluid_vol_mm3(k))), PB.peak_dV_mV(k)*1e3, 'uV', ...
        'noiseless; baseline = healthy head on the same montage');
    R = addrow(R, 'b: peak |dV| relative to that reading', ...
        sprintf('primary | pooled 3 montages | burden %s', bfmt(b)), ...
        PB.peak_dV_rel(k)*1e6, 'ppm', ...
        'peak |dV| / |V0| on the channel where the peak occurs, parts per million');
    for m = MONTAGES
        km = k & PB.montage == m;
        R = addrow(R, 'b: peak |dV| (fluid minus healthy)', ...
            sprintf('primary | %s | burden %s', m, bfmt(b)), PB.peak_dV_mV(km)*1e3, 'uV', ...
            'noiseless; baseline = healthy head on the same montage');
    end
end

%% (c) Effect of each axis (dB, paired)
fprintf('\n=============== (c) PAIRED AXIS EFFECTS (dB) ===============\n');

% Fluid volume
[G, ~, ~, ~] = findgroups(PB.montage, PB.anchor, PB.depth_mm);
fx = nan(max(G), 1); fr = fx;
for g = 1:max(G)
    k = find(G == g); b = PB.burden_frac(k); p = PB.peak_dV_mV(k); v = PB.fluid_vol_mm3(k);
    fx(g) = 20*log10(p(b == max(b)) / p(b == min(b)));
    fr(g) = v(b == max(b)) / v(b == min(b));
end
R = addrow(R, 'c: effect of fluid volume', 'primary | 10% vs 0.25% burden', fx, 'dB', ...
    sprintf(['pair: same montage, anchor, depth; 20log10(peak at 10%% / peak at 0.25%%). ' ...
    'Fluid-volume ratio %.1f-%.1f (%.2f dB for a 1:1 volume law)'], min(fr), max(fr), 20*log10(40)));

% Inclusion shape
[~, loc] = ismember(replace(PB.key, PRIMARY, ALTRULE), PBL.key);
sx = 20*log10(PBL.peak_dV_mV(loc) ./ PB.peak_dV_mV);
R = addrow(R, 'c: effect of inclusion shape', 'primary geometry (scenario B) | L=0 vs L=1/3', ...
    sx, 'dB', ['pair: identical configuration, only the mixing rule differs; 20log10(peak L=0 / peak L=1/3). ' ...
    'The spread is almost entirely the burden grid (see by-burden rows). Element-level conductivity-increment ' ...
    'ratio from the mixing formula: 16.4 dB as f->0, 15.7 dB at f = 10% (a2_mix_conductivity_L:75-79 ' ...
    'quotes the 0-10% chord slopes as initial slopes)']);
for b = BURDENS
    k = abs(PB.burden_frac - b) < 1e-12;
    R = addrow(R, 'c: effect of inclusion shape, by burden', ...
        sprintf('primary geometry (scenario B) | L=0 vs L=1/3 | burden %s', bfmt(b)), sx(k), 'dB', ...
        'same twin pairing, one burden level; the effect shrinks as burden rises because L=1/3 mixing is convex in f');
end
AM = I(I.scenario == "A" & I.rule == PRIMARY, :); AL = I(I.scenario == "A" & I.rule == ALTRULE, :);
[~, loc] = ismember(replace(AM.key, PRIMARY, ALTRULE), AL.key);
R = addrow(R, 'c: effect of inclusion shape', 'single-space bound (scenario A) | L=0 vs L=1/3', ...
    20*log10(AL.peak_dV_mV(loc) ./ AM.peak_dV_mV), 'dB', ...
    'scenario A only; a bound, not pooled with scenario B');

% Depth
[G, ~, ga, ~] = findgroups(PB.montage, PB.anchor, PB.burden_frac);   % ga = anchor of each group
dx = nan(max(G), 1); dsp = dx; dper = dx; mono = false(max(G), 1); span = dx;
for g = 1:max(G)
    k = find(G == g); [d, ix] = sort(PB.depth_mm(k)); p = PB.peak_dV_mV(k(ix));
    dx(g)   = 20*log10(p(1) / p(end));
    dsp(g)  = 20*log10(max(p) / min(p));
    span(g) = d(end) - d(1);
    dper(g) = dx(g) / span(g) * 10;
    mono(g) = all(diff(p) < 0);
end
R = addrow(R, 'c: effect of depth', 'primary | shallowest (25 mm) vs deepest reachable', dx, 'dB', ...
    sprintf(['pair: same anchor, burden, montage; 20log10(peak at 25 mm / peak at deepest). ' ...
    'Span 10-30 mm depending on anchor. Signal falls monotonically with depth in %d of %d groups'], ...
    sum(mono), numel(mono)));
R = addrow(R, 'c: effect of depth', 'primary | max/min spread over reachable depths', dsp, 'dB', ...
    'same groups; unsigned spread, the definition used in the Methods');
R = addrow(R, 'c: effect of depth, per 10 mm', 'primary | shallowest vs deepest, normalised', ...
    dper, 'dB per 10 mm', 'effect divided by the depth span of each group, x10');
for a = ["Cz","Fz","C3","Pz","T7","O1"]
    k = ga == a;
    R = addrow(R, 'c: effect of depth, by anchor', ...
        sprintf('primary | %s | 25 mm vs %d mm', a, 25 + median(span(k))), dx(k), 'dB', ...
        'same pairing, one anchor; shows the depth effect depends on how deep the line can go');
end
Cc = PB(PB.depth_mm <= CORE_MAX_DEPTH, :);
[G, ~, ~, ~] = findgroups(Cc.montage, Cc.anchor, Cc.burden_frac);
dc = nan(max(G), 1);
for g = 1:max(G)
    k = find(G == g); p25 = Cc.peak_dV_mV(k(Cc.depth_mm(k) == 25)); p35 = Cc.peak_dV_mV(k(Cc.depth_mm(k) == 35));
    dc(g) = 20*log10(p25 / p35);
end
R = addrow(R, 'c: effect of depth', 'primary | 25 vs 35 mm (balanced core, all 6 anchors)', dc, 'dB', ...
    'pair: same anchor, burden, montage; the only depth span every anchor can reach');

% Position
[G, ~, ~, ~] = findgroups(Cc.montage, Cc.depth_mm, Cc.burden_frac);
px = nan(max(G), 1); na = px; hi = strings(max(G), 1); lo = hi;
for g = 1:max(G)
    k = find(G == g); p = Cc.peak_dV_mV(k);
    px(g) = 20*log10(max(p) / min(p)); na(g) = numel(k);
    [~, i1] = max(p); [~, i2] = min(p); hi(g) = Cc.anchor(k(i1)); lo(g) = Cc.anchor(k(i2));
end
assert(all(na == 6), 'Balanced core does not have all six anchors in every group.');
R = addrow(R, 'c: effect of position', 'primary | spread across 6 anchors, 25-35 mm', px, 'dB', ...
    sprintf(['pair: same depth, burden, montage; 20log10(max/min over the 6 anchors). ' ...
    'Strongest anchor: %s; weakest: %s'], tally(hi), tally(lo)));
[G, ~, ~, ~] = findgroups(PB.montage, PB.depth_mm, PB.burden_frac);
pa = nan(max(G), 1); naa = pa;
for g = 1:max(G)
    p = PB.peak_dV_mV(G == g); naa(g) = numel(p); pa(g) = 20*log10(max(p) / min(p));
end
R = addrow(R, 'c: effect of position', 'primary | spread across reachable anchors, all depths', ...
    pa(naa >= 2), 'dB', sprintf(['secondary: deeper groups hold only %d-%d anchors, so this ' ...
    'mixes spread with how many anchors exist; the balanced-core row is the clean one'], ...
    min(naa(naa >= 2)), max(naa)));

% Region size at fixed fluid volume
LV = [0 1/3 1/2]; LN = ["L=0", "L=1/3", "L=1/2"];
for q = 1:3
    S = Z(abs(Z.L - LV(q)) < 1e-9, :);
    [G, ~, ~, ~] = findgroups(S.anchor, S.depth_target_mm, S.fluid_vol_target_mm3);
    zx = nan(max(G), 1); zs = zx;
    for g = 1:max(G)
        k = find(G == g); [~, ix] = sort(S.diam_mm(k)); p = S.peak_dV_mV(k(ix));
        zx(g) = 20*log10(max(p) / min(p)); zs(g) = 20*log10(p(end) / p(1));
    end
    R = addrow(R, 'c/d: effect of region size at fixed fluid volume', ...
        sprintf('size sweep | %s | spread over 12-32 mm', LN(q)), zx, 'dB', ...
        sprintf(['pair: same anchor (Cz, Fz), depth (30, 40 mm), fluid volume (100/200/400 mm3); ' ...
        '20log10(max/min over 5 diameters). Signed 32 mm vs 12 mm: median %+.2f dB ' ...
        '(negative = spreading the fluid wider weakens the signal)'], median(zs)));
end

% Montage
for is = 1:numel(SNR_DB)
    Ps = P(P.snr_db == SNR_DB(is), :);
    [G, ~, ~, ~] = findgroups(Ps.anchor, Ps.depth_mm, Ps.burden_frac);
    mx = nan(max(G), 1); nm3 = mx; win = strings(max(G), 1);
    for g = 1:max(G)
        k = find(G == g); mg = Ps.margin_2frame_p95_dB(k); nm3(g) = numel(k);
        mx(g) = max(mg) - min(mg); [~, i1] = max(mg); win(g) = Ps.montage(k(i1));
    end
    assert(all(nm3 == 3), 'Montage groups are not complete.');
    R = addrow(R, 'c: effect of montage (on margin)', ...
        sprintf('primary | spread across 3 montages | %d dB', SNR_DB(is)), mx, 'dB', ...
        sprintf(['pair: identical object (anchor, depth, burden); max-min margin over the 3 ' ...
        'montages. Best montage: %s'], tally(win)));
end
[G, ~, ~, ~] = findgroups(PB.anchor, PB.depth_mm, PB.burden_frac);
mr = nan(max(G), 1); winr = strings(max(G), 1);
for g = 1:max(G)
    k = find(G == g); p = PB.peak_dV_mV(k);
    mr(g) = 20*log10(max(p) / min(p)); [~, i1] = max(p); winr(g) = PB.montage(k(i1));
end
R = addrow(R, 'c: effect of montage (raw peak, before noise)', ...
    'primary | spread across 3 montages', mr, 'dB', ...
    sprintf(['same pairing on noiseless peak |dV|. Largest raw peak: %s. The gap between this ' ...
    'and the margin rows is the frame-norm noise scaling (a2_noise_model:80-90), an open decision'], ...
    tally(winr)));

% Monte Carlo stability of the p95 comparator, which explains why the montage row changes with SNR
for m = MONTAGES
    k = F.montage == m;
    x95 = 20*log10(F.peak_noise_p95_2frame_mV(k) .* F.snr_linear(k));
    x50 = 20*log10(F.peak_noise_only_2frame_mV(k) .* F.snr_linear(k));
    R = addval(R, 'noise comparator: p95 wander across SNR levels', char("noise only | " + m), ...
        sum(F.n_real(k)), max(x95) - min(x95), 'dB', sprintf(['single value in median column. ' ...
        'SNR-normalised p95 noise-only peak should be identical at every SNR; it varies by this ' ...
        'much = Monte Carlo error of a p95 from 200 realisations. Median comparator varies %.2f dB'], ...
        max(x50) - min(x50)));
end
for j = 1:height(F)
    R = addval(R, 'noise comparator: noise-only peak, two-frame', ...
        sprintf('noise only | %s | %d dB', F.montage(j), F.snr_db(j)), F.n_real(j), ...
        F.peak_noise_only_2frame_mV(j)*1e3, 'uV', sprintf(['median over %d realisations in the ' ...
        'median column. p95 = %.4g uV (the detection threshold). Per-channel sigma %.4g uV. SD of the ' ...
        'noise-only peak is not in the saved outputs'], F.n_real(j), ...
        F.peak_noise_p95_2frame_mV(j)*1e3, F.sigma_2frame_mV(j)*1e3));
end

%% (e) Regression of log10 peak |dV| on log10 fluid volume
fprintf('\n=============== (e) FLUID-VOLUME REGRESSION ===============\n');
for m = [MONTAGES, "pooled"]
    if m == "pooled", k = true(height(PB), 1); else, k = PB.montage == m; end
    x = log10(PB.fluid_vol_mm3(k)); y = log10(PB.peak_dV_mV(k)); n = numel(x);
    X = [ones(n, 1) x]; bcoef = X \ y; res = y - X*bcoef;
    r2 = 1 - sum(res.^2) / sum((y - mean(y)).^2);
    se = sqrt(sum(res.^2) / (n - 2) / sum((x - mean(x)).^2));
    R = addval(R, 'e: log-log slope, peak |dV| vs fluid volume', char("primary | " + m), n, ...
        bcoef(2), '-', sprintf(['OLS point estimate in median column. SE %.4f; approx 95%% CI ' ...
        '%.3f-%.3f (+/-1.96 SE). Slope 1 = peak directly proportional to fluid volume'], ...
        se, bcoef(2) - 1.96*se, bcoef(2) + 1.96*se));
    R = addval(R, 'e: R^2, variance explained by fluid volume', char("primary | " + m), n, ...
        r2, 'fraction', sprintf(['single value in median column; r = %.4f. Pooled row also ' ...
        'absorbs the montage offsets, so it is lower than the per-montage rows'], sqrt(r2)));
end

%% (f) Single-space bound (scenario A, L = 1/3)
fprintf('\n=============== (f) SINGLE-SPACE BOUND (scenario A) ===============\n');
A = N(N.scenario == "A", :);
for s = [45 48]
    k = A.rule == PRIMARY & A.snr_db == s;
    R = addrow(R, 'f: single-space bound margin', ...
        sprintf('scenario A (3 and 5 mm) | L=1/3 | 3 montages | %d dB', s), ...
        A.margin_2frame_p95_dB(k), 'dB', ['BOUND ONLY: under-resolved inclusion over-responds, so ' ...
        '"undetectable here" is safe and "detectable here" would not be. max = best case: ' ...
        bestcase(A(k,:))]);
end
k = A.rule == ALTRULE & A.snr_db == 48;
R = addrow(R, 'f: single-space bound margin, RETIRED settings', ...
    'scenario A | L=0 | two-frame | MEDIAN comparator | 48 dB', A.margin_2frame_dB(k), 'dB', ...
    ['provenance only: reproduces the old "17 dB below" figure (max). Not a primary number. max = ' ...
    bestcase(A(k,:), 'margin_2frame_dB')]);

%% Example configuration for the healthy versus fluid figure
ke = PB.montage == EX.montage & PB.anchor == EX.anchor & PB.depth_mm == EX.depth_mm & ...
     abs(PB.burden_frac - EX.burden_frac) < 1e-12;
assert(sum(ke) == 1, 'Example configuration not found exactly once.');
ex  = PB(ke, :);
v0  = VB.(EX.montage);
dv  = DVM.(EX.montage)(char(ex.key));
Pm  = P(P.montage == EX.montage, :);
exm = @(s) Pm.margin_2frame_p95_dB(Pm.key == ex.key & Pm.snr_db == s);
pr48 = 100 * mean(Pm.margin_2frame_p95_dB(Pm.snr_db == 48) <= exm(48));
[pkv, ipk] = max(abs(dv));
R = addval(R, 'fig example: peak |dV|', char(ex.key), numel(dv), pkv*1e6, 'uV', ...
    sprintf(['single value. Channel %d; healthy reading there %.1f uV (%.0f ppm change). ' ...
    'Fluid %.1f mm3, depth achieved %.2f mm'], ipk, abs(v0(ipk))*1e6, pkv/abs(v0(ipk))*1e6, ...
    ex.fluid_vol_mm3, ex.depth_actual_mm));
for s = SNR_DB
    R = addval(R, 'fig example: margin above noise alone', sprintf('%s | %d dB', ex.key, s), 1, ...
        exm(s), 'dB', 'single value');
end
R = addval(R, 'fig example: percentile of its 48 dB margin', ...
    char("within primary | " + EX.montage), sum(Pm.snr_db == 48), pr48, '%', ...
    'single value; 50 = exactly typical for this montage');

% Folder for the Results figures drawn by later scripts
if ~exist(figd, 'dir'), mkdir(figd); end

%% Save and print the summary table
numcols = {'n','median','sd','iqr_lo','iqr_hi','min','max'};
assert(all(varfun(@(c) isa(c, 'double'), R(:, numcols), 'OutputFormat', 'uniform')), ...
    'A numeric column is not double; values would be rounded on write.');
writetable(R, outcsv);
% Check that the saved file matches what was computed.
Rb = readtable(outcsv, 'TextType', 'string');
e = max(abs(Rb{:, numcols} - R{:, numcols}), [], 'all', 'omitnan');
assert(height(Rb) == height(R) && e < 1e-9 * max(1, max(abs(R{:, numcols}), [], 'all', 'omitnan')), ...
    'CSV read-back does not match the computed table.');
fprintf('\nCSV read-back: %d rows, max abs difference %.2e (PASS)\n', height(Rb), e);
fprintf('\n=============== HEADLINE NUMBERS: median +/- SD (n) ===============\n');
hsel = ~startsWith(R.quantity, "check") & ~startsWith(R.quantity, "noise comparator: noise-only") ...
     & ~contains(R.quantity, "by burden") & ~contains(R.quantity, "by anchor") ...
     & ~contains(R.quantity, "per reading") ...
     & ~(startsWith(R.quantity, "b: peak") & ~contains(R.subset, "pooled"));
H = R(hsel, :);
for j = 1:height(H)
    fprintf('  %-48s %-62s %9.3f +/- %-8.3f (n=%d) %s\n', H.quantity(j), H.subset(j), ...
        H.median(j), H.sd(j), H.n(j), H.units(j));
end
fprintf('\nSAVED %s   (%d rows)\n', outcsv, height(R));

%% Local functions
function R = addrow(R, quantity, subset, x, units, note)
% Summarise a sample, append it as one row and print it.
x = double(x(:));
[g1, bc] = shape_stats(x);
flag = '';
if ~isnan(g1) && abs(g1) > 1, flag = [flag ' SKEWED (|G1|>1): the IQR is the better spread.']; end
if ~isnan(bc) && bc > 5/9,     flag = [flag ' FLAT OR MULTI-PEAKED (BC>5/9): SD spans clusters, not scatter about one value.']; end
note = sprintf('%s | skew G1 %.2f, bimodality BC %.2f.%s', note, g1, bc, flag);
row = table(string(quantity), string(subset), numel(x), median(x), std(x), pctl(x, 25), ...
    pctl(x, 75), min(x), max(x), string(units), string(note), 'VariableNames', ...
    {'quantity','subset','n','median','sd','iqr_lo','iqr_hi','min','max','units','note'});
R = [R; row];
fprintf('  %-44s | %-58s n=%4d  med %9.3f  sd %7.3f  IQR [%9.3f, %9.3f]  range [%9.3f, %9.3f] %s%s\n', ...
    quantity, subset, numel(x), median(x), std(x), pctl(x, 25), pctl(x, 75), min(x), max(x), ...
    units, flag);
end

function R = addval(R, quantity, subset, n, value, units, note)
% Append a single-value row (value stored in the median column).
n = double(n); value = double(value);
row = table(string(quantity), string(subset), n, value, NaN, NaN, NaN, NaN, NaN, ...
    string(units), string(note), 'VariableNames', ...
    {'quantity','subset','n','median','sd','iqr_lo','iqr_hi','min','max','units','note'});
R = [R; row];
fprintf('  %-44s | %-58s n=%4d  value %10.4f %s\n', quantity, subset, n, value, units);
end

function [R, allok] = addcheck(R, allok, name, value, ok, detail)
% Record an integrity check in the output table and on screen.
if ok, v = 'PASS'; else, v = '**FAIL**'; end
% Convert to double, so an integer value does not change the type of the whole column.
value = double(value);
row = table("check", string(name), NaN, value, NaN, NaN, NaN, NaN, NaN, "-", ...
    string(sprintf('%s: %s', v, detail)), 'VariableNames', ...
    {'quantity','subset','n','median','sd','iqr_lo','iqr_hi','min','max','units','note'});
R = [R; row];
allok = allok && ok;
fprintf('  %-8s %-62s value %-12.4g %s\n', v, name, value, detail);
end

function s = bestcase(T, col)
% Describe the configuration with the largest value of T.(col).
if nargin < 2, col = 'margin_2frame_p95_dB'; end
[~, i] = max(T.(col));
s = sprintf('%s, %s, %g mm, size %g mm, burden %g%%', T.montage(i), T.anchor(i), ...
    T.depth_mm(i), T.size_mm(i), T.burden_frac(i)*100);
end

function s = tally(x)
% Counts of each value in a string vector, most frequent first.
[u, ~, j] = unique(x); c = accumarray(j, 1); [c, o] = sort(c, 'descend'); u = u(o);
s = strjoin(compose('%s %d', u, c), ', ');
end

function q = pctl(x, p)
% Percentile without the Statistics Toolbox, using the same interpolation as
% MATLAB's prctile.
x = sort(x(:)); n = numel(x);
assert(n > 0, 'pctl:empty', 'Empty sample: a subset selected no rows.');
if n == 1, q = x; return, end
pos = p/100*n + 0.5;
if pos <= 1, q = x(1); return, end
if pos >= n, q = x(n); return, end
lo = floor(pos); q = x(lo) + (pos - lo)*(x(lo + 1) - x(lo));
end

function [G1, BC] = shape_stats(x)
% Bias-corrected skewness G1 and Sarle's bimodality coefficient.
%   BC = (G1^2 + 1) / (G2 + 3(n-1)^2/((n-2)(n-3))), G2 = bias-corrected excess
%   kurtosis. A uniform distribution gives BC = 5/9.
n = numel(x);
if n < 4 || std(x) == 0, G1 = NaN; BC = NaN; return, end
d = x - mean(x); m2 = mean(d.^2);
g1 = mean(d.^3) / m2^1.5; g2 = mean(d.^4) / m2^2 - 3;
G1 = sqrt(n*(n - 1)) / (n - 2) * g1;
G2 = (n - 1) / ((n - 2)*(n - 3)) * ((n + 1)*g2 + 6);
BC = (G1^2 + 1) / (G2 + 3*(n - 1)^2 / ((n - 2)*(n - 3)));
end

function r = thesis_root()
% Path to the repository root.
r = fileparts(mfilename('fullpath'));
while ~exist(fullfile(r, 'README.md'), 'file')
    p = fileparts(r);
    if isempty(p) || strcmp(p, r)
        error('thesis_root:noRoot', 'Could not find the repository root (README.md).');
    end
    r = p;
end
end
