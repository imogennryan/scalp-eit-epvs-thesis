%A2_ADD_NOISE  Add Gaussian noise to the simulated dataset and measure the peak against noise alone.

% Aim served: Aim 2.
%
% Two noise cases are computed. Two-frame (main case): the baseline and the
% perturbed frame both carry noise, so the noise on the difference is sqrt(2)
% larger. One-frame: the baseline is treated as exactly known.
% One noise ensemble per montage and SNR is reused for every configuration,
% so configurations are compared on the same noise.
% A peak over many channels is large even with no signal, so each peak is
% compared with the peak that noise alone gives (median and 95th percentile).
%
% Inputs : data/aim2/a2_sweep_*.mat, results/aim2/a2_sweep_index.csv,
%          a2_noise_model
% Outputs: results/aim2/a2_noise_summary.csv   one row per configuration x SNR
%          results/aim2/a2_noise_floor.csv     one row per montage x SNR
% Imogen Ryan - Honours thesis.

clear; close all force;
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(fileparts(here),'aim1_head_model'));

outd = thesis_dir('results','aim2');
datd = thesis_dir('data','aim2');

SNR_DB   = [45 48 77 95];    % instrument noise levels (dB), Goren 2018, McDermott 2020, Ouypornkochagorn 2022
N_REAL   = 200;              % realisations per (montage, SNR)
SEED_A   = 0;                % baseline frame noise
SEED_B   = 1;                % second frame noise, for the two-frame case

fprintf('\n############### A2: ADD NOISE TO THE DATASET ###############\n');

idx = readtable([outd 'a2_sweep_index.csv'], 'TextType','string');
fprintf('index: %d configurations\n', height(idx));

montages = unique(idx.montage, 'stable');
rows  = {};
floorrows = {};

for im = 1:numel(montages)
    mtg = montages(im);
    S  = load([datd sprintf('a2_sweep_%s.mat', mtg)], 'DV', 'vb');
    vb = S.vb;  nm = numel(vb);

    sel = find(idx.montage == mtg);
    fprintf('\n--- %s : %d measurements, %d configurations ---\n', mtg, nm, numel(sel));

    % Two independent noise ensembles, nm x N_REAL x numel(SNR_DB).
    [NA, ninfo] = a2_noise_model(vb, SNR_DB, N_REAL, SEED_A);
    NB          = a2_noise_model(vb, SNR_DB, N_REAL, SEED_B);

    for is = 1:numel(SNR_DB)
        n1 = NA(:,:,is);            % one-frame  : noise on the difference
        n2 = NB(:,:,is) - n1;       % two-frame  : (frame2 - frame1)

        % Noise standard deviation measured from the ensemble.
        sig1 = std(n1(:));
        sig2 = std(n2(:));

        % Peak of noise alone over nm channels. The 95th percentile is the
        % detection threshold (5% false-positive rate per frame); the median
        % is the typical noise peak.
        pkd1 = max(abs(n1), [], 1);   pkd2 = max(abs(n2), [], 1);
        pk_noise1 = median(pkd1);     pk_noise2 = median(pkd2);
        pk95_1    = prctile_local(pkd1, 95);
        pk95_2    = prctile_local(pkd2, 95);

        floorrows{end+1} = table(mtg, nm, SNR_DB(is), ninfo.snr_linear(is), ...
            sig1*1e3, sig2*1e3, pk_noise1*1e3, pk_noise2*1e3, ...
            pk95_1*1e3, pk95_2*1e3, 20*log10(pk95_2/pk_noise2), ...
            sig2/sig1, 20*log10(sig2/sig1), N_REAL, ...
            'VariableNames', {'montage','n_meas','snr_db','snr_linear', ...
            'sigma_1frame_mV','sigma_2frame_mV','peak_noise_only_1frame_mV', ...
            'peak_noise_only_2frame_mV','peak_noise_p95_1frame_mV', ...
            'peak_noise_p95_2frame_mV','p95_over_median_dB', ...
            'two_over_one_frame','penalty_dB','n_real'}); %#ok<SAGROW>

        for jj = 1:numel(sel)
            r  = sel(jj);
            ky = char(idx.key(r));
            dv = S.DV(ky);

            pk0 = max(abs(dv));                       % noiseless peak

            pk1 = median(max(abs(dv + n1), [], 1));   % one-frame noisy peak
            pk2 = median(max(abs(dv + n2), [], 1));   % two-frame noisy peak

            rows{end+1} = {ky, mtg, idx.scenario(r), idx.anchor(r), ...
                idx.depth_mm(r), idx.size_mm(r), idx.burden_frac(r), ...
                idx.rule(r), SNR_DB(is), ...
                pk0*1e3, pk1*1e3, pk2*1e3, ...
                pk_noise1*1e3, pk_noise2*1e3, ...
                sig1*1e3, sig2*1e3, ...
                20*log10(pk0/pk_noise1), 20*log10(pk0/pk_noise2), ...
                20*log10(pk0/pk95_1), 20*log10(pk0/pk95_2), ...
                pk0/sig1, pk0/sig2, ...
                norm(dv)/(norm(vb)/ninfo.snr_linear(is))}; %#ok<SAGROW>
        end
        fprintf('  %2d dB : sigma %8.2e mV (1-frame) %8.2e mV (2-frame), ', ...
            SNR_DB(is), sig1*1e3, sig2*1e3);
        fprintf('noise-only peak %8.2e mV\n', pk_noise2*1e3);
    end
end

%% Assemble and save
V = {'key','montage','scenario','anchor','depth_mm','size_mm','burden_frac', ...
     'rule','snr_db','peak_dV_mV','peak_noisy_1frame_mV','peak_noisy_2frame_mV', ...
     'peak_noise_only_1frame_mV','peak_noise_only_2frame_mV', ...
     'sigma_1frame_mV','sigma_2frame_mV', ...
     'margin_1frame_dB','margin_2frame_dB', ...
     'margin_1frame_p95_dB','margin_2frame_p95_dB', ...
     'peak_over_sigma_1frame','peak_over_sigma_2frame','dnr_norm'};
T = cell2table(vertcat(rows{:}), 'VariableNames', V);
writetable(T, [outd 'a2_noise_summary.csv']);

F = vertcat(floorrows{:});
writetable(F, [outd 'a2_noise_floor.csv']);

%% Report
% Margin is the peak signal against the noise-only peak, in dB. Positive means
% the peak is above it, negative is the shortfall.
fprintf('\n=================== MARGIN SUMMARY ===================\n');
fprintf('Peak signal vs the p95 noise-only peak (5%% false-positive rate per\n');
fprintf('frame), two-frame. Negative = below. L=1/3 is the physically\n');
fprintf('defensible shape case; L=0 (linear) is the aligned-tube extreme.\n\n');
RL = ["maxwell_garnett","linear"];
RN = ["L=1/3 (spheres/random)","L=0   (aligned tubes)"];
fprintf('  %-6s %-24s %10s %10s %10s\n','SNR','shape case','best','median','worst');
for is = 1:numel(SNR_DB)
    for q = 1:2
        m = T.snr_db == SNR_DB(is) & T.rule == RL(q);
        fprintf('  %-6d %-24s %10.1f %10.1f %10.1f\n', SNR_DB(is), RN(q), ...
            max(T.margin_2frame_p95_dB(m)), median(T.margin_2frame_p95_dB(m)), ...
            min(T.margin_2frame_p95_dB(m)));
    end
end

% Margins by burden and depth. The grid is unbalanced because unreachable
% depths are left out, so slices are reported instead of overall fractions.
fprintf('\n============ CONDITIONAL MARGINS (L=1/3, two-frame, p95) ============\n');
fprintf('Scenario B, 19-electrode montage. Margin in dB by burden and depth.\n\n');
sel0 = T.rule=="maxwell_garnett" & T.scenario=="B" & T.montage=="m19_1020";
for is = [2 4]                                    % 48 dB and 95 dB
    fprintf('  --- %d dB ---\n', SNR_DB(is));
    fprintf('  %-10s', 'burden');
    for d = [25 35 45 55], fprintf('%9s', sprintf('%dmm', d)); end
    fprintf('\n');
    for bf = [0.0025 0.01 0.05 0.10]
        fprintf('  %-10s', sprintf('%.2f%%', bf*100));
        for d = [25 35 45 55]
            k = sel0 & T.snr_db==SNR_DB(is) & abs(T.burden_frac-bf)<1e-9 & T.depth_mm==d;
            if any(k), fprintf('%9.1f', median(T.margin_2frame_p95_dB(k)));
            else,      fprintf('%9s', '-'); end
        end
        fprintf('\n');
    end
    fprintf('\n');
end

fprintf('For reference, the same cells under L=0 (aligned tubes, +15.7 dB on\n');
fprintf('the mixing slope) sit about 14 dB higher.\n');

fprintf('\nOne-frame minus two-frame margin (should be ~%.2f dB throughout): %.2f to %.2f dB\n', ...
    20*log10(sqrt(2)), min(T.margin_1frame_dB - T.margin_2frame_dB), ...
    max(T.margin_1frame_dB - T.margin_2frame_dB));
fprintf('p95 comparator is %.1f to %.1f dB stricter than the median comparator.\n', ...
    min(T.margin_2frame_dB - T.margin_2frame_p95_dB), ...
    max(T.margin_2frame_dB - T.margin_2frame_p95_dB));

fprintf('\nSAVED %s   (%d rows)\n', [outd 'a2_noise_summary.csv'], height(T));
fprintf('SAVED %s   (%d rows)\n', [outd 'a2_noise_floor.csv'], height(F));

function q = prctile_local(x, p)
% Percentile without the Statistics Toolbox, using the same interpolation as
% MATLAB's prctile.
x = sort(x(:)); n = numel(x);
if n == 1, q = x; return, end
pos = p/100*n + 0.5;
if pos <= 1, q = x(1); return, end
if pos >= n, q = x(n); return, end
lo = floor(pos); q = x(lo) + (pos-lo)*(x(lo+1)-x(lo));
end

function d = thesis_dir(varargin)
% Path to a folder inside this repository.
r = fileparts(mfilename('fullpath'));
while ~exist(fullfile(r,'README.md'),'file')
    p = fileparts(r);
    if isempty(p) || strcmp(p,r)
        error('thesis_dir:noRoot','Could not find the repository root (README.md).');
    end
    r = p;
end
d = [fullfile(r, varargin{:}) filesep];
end
