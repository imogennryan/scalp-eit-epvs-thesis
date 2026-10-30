%A1_GENERATE_BASELINE  Solve and save the healthy baseline for each montage.

% Aim served: Aim 1. Every Aim 2 perturbation is compared against this baseline.
%
% The baseline is the four-layer model with no perturbation. This is a modelling
% reference, not a real healthy brain, which usually has a few small
% perivascular spaces. Each of the three montages gets its own baseline.
%
% Inputs : a1_build_head_model, models/montages/a1_montages.mat
% Outputs: models/baseline/a1_baseline_<montage>.mat (model and provenance),
%          data/aim1/a1_baseline_meas_<montage>.mat (measurements),
%          results/aim1/a1_baseline_summary.csv
% Imogen Ryan - Honours thesis.

clear; close all force;
addpath(fileparts(mfilename('fullpath')));
eidors_msg('log_level',1);
rng('default');                       % reproducibility

mdld = thesis_dir('models','baseline'); if ~exist(mdld,'dir'), mkdir(mdld); end
datd = thesis_dir('data','aim1');       if ~exist(datd,'dir'), mkdir(datd); end
outd = thesis_dir('results','aim1');    if ~exist(outd,'dir'), mkdir(outd); end

montages = {'m19_1020','m32_1010','m16_twoplane'};
T = table();

fprintf('\n############ A1: HEALTHY BASELINE GENERATION ############\n');
for k = 1:numel(montages)
    t0 = tic;
    [img, info] = a1_build_head_model(struct('montage', montages{k}));
    ne = numel(img.fwd_model.electrode);

    vh = fwd_solve(img);
    v  = vh.meas;
    t_solve = toc(t0);

    % Checks that can fail, saved with the model
    chk = struct();
    chk.n_meas_expected = ne*(ne-3);
    chk.n_meas          = numel(v);
    chk.count_ok        = chk.n_meas == chk.n_meas_expected;
    chk.all_finite      = all(isfinite(v));
    chk.no_exact_zeros  = ~any(v==0);
    % Every electrode is a complete electrode model, with no orphan nodes
    B = img.fwd_model.boundary; npoint = 0; norphan = 0;
    for e = 1:ne
        en = img.fwd_model.electrode(e).nodes;
        bi = find_electrode_bdy(B, img.fwd_model.nodes, en);
        if isempty(bi), npoint = npoint + 1; end
        norphan = norphan + numel(setdiff(en, reshape(B(bi,:),1,[])));
    end
    chk.all_CEM   = (npoint==0) && (norphan==0);
    % Reciprocity: swapping the drive pair and the measurement pair must give
    % the same reading, so the full ne x ne matrix M must be symmetric.
    % This shows the equations were solved consistently, not that the model is correct.
    [s2,m2] = mk_stim_patterns(ne,1,'{ad}','{ad}',{'meas_current'},info.opt.amp);
    i2 = img; i2.fwd_model.stimulation = s2; i2.fwd_model.meas_select = m2;
    M  = reshape(fwd_solve(i2).meas, ne, ne);
    chk.reciprocity = max(max(abs(M-M.')))/max(abs(M(:)));
    chk.recip_ok    = chk.reciprocity < 1e-10;
    chk.ALL_PASS = chk.count_ok && chk.all_finite && chk.no_exact_zeros && ...
                   chk.all_CEM && chk.recip_ok;

    % Provenance, saved inside the .mat with the model
    prov = struct();
    prov.generated       = datestr(now,'yyyy-mm-dd HH:MM:SS');
    prov.script          = mfilename('fullpath');
    prov.mesh            = info.opt.mesh;
    prov.n_elements      = size(img.fwd_model.elems,1);
    prov.montage         = montages{k};
    prov.n_electrodes    = ne;
    prov.elec_area_mm2   = info.elec_area_m2*1e6;
    prov.elec_diam_mm    = 2*sqrt(median(info.elec_area_m2)/pi)*1e3;
    prov.elec_snap_mm    = info.elec_snap_err_mm;
    prov.z_contact       = info.opt.z_contact;
    prov.z_contact_src   = 'McDermott 2020: 1 kOhm on 10 mm electrodes';
    prov.current_A       = info.opt.amp;
    prov.current_src     = '100 uA rms safety limit up to 1 kHz, McDermott 2020 citing ISO 2002';
    prov.sigma           = info.sigma;
    prov.sigma_src       = info.sigma_source;
    prov.units           = info.units;
    prov.stim_protocol   = 'adjacent drive, adjacent measure, no_meas_current';
    prov.checks          = chk;
    prov.drive_frequency = ['100 Hz. Top of the band where ' ...
                            'every conductivity sits inside its own source ' ...
                            '(McDermott 5-100 Hz; Baumann 10 Hz-10 kHz). ' ...
                            '100 uA limit applies below 1 kHz.'];
    prov.known_mismatch  = ['Research-system noise floors (Ouypornkochagorn ' ...
                            '2023, 77-95 dB) are measured at 10 kHz. Comparing ' ...
                            'against them is a cross-frequency comparison.'];

    save([mdld sprintf('a1_baseline_%s.mat',montages{k})], 'img','info','prov','-v7.3');
    save([datd sprintf('a1_baseline_meas_%s.mat',montages{k})], 'v','prov');

    fprintf('\n--- %s : %d electrodes ---\n', montages{k}, ne);
    fprintf('  measurements   : %d (expected %d)   %s\n', chk.n_meas, chk.n_meas_expected, pf(chk.count_ok));
    fprintf('  all finite     : %s   no exact zeros: %s\n', pf(chk.all_finite), pf(chk.no_exact_zeros));
    fprintf('  all CEM        : %s  (%d point electrodes, %d orphan nodes)\n', pf(chk.all_CEM), npoint, norphan);
    fprintf('  reciprocity    : %.3e   %s\n', chk.reciprocity, pf(chk.recip_ok));
    fprintf('  voltages       : RMS %.4f mV | max %.4f mV | min %.4f mV\n', ...
        rms(v)*1e3, max(abs(v))*1e3, min(abs(v))*1e3);
    fprintf('  transfer Z     : median %.2f Ohm | max %.2f Ohm\n', ...
        median(abs(v))/info.opt.amp, max(abs(v))/info.opt.amp);
    fprintf('  solve time     : %.2f s\n', t_solve);
    fprintf('  OVERALL        : %s\n', pf(chk.ALL_PASS));

    T = [T; table(string(montages{k}), ne, chk.n_meas, chk.reciprocity, ...
        rms(v)*1e3, max(abs(v))*1e3, median(abs(v))/info.opt.amp, ...
        prov.elec_diam_mm, max(info.elec_snap_err_mm), chk.ALL_PASS, ...
        'VariableNames',{'montage','n_elec','n_meas','reciprocity','rms_mV', ...
        'max_mV','medZ_Ohm','elec_diam_mm','snap_max_mm','all_checks_pass'})];
end

writetable(T,[outd 'a1_baseline_summary.csv']);
fprintf('\n############ BASELINES SAVED ############\n');
fprintf('models : %s\n', mdld);
fprintf('data   : %s\n', datd);
fprintf('summary: %sa1_baseline_summary.csv\n', outd);
if all(T.all_checks_pass)
    fprintf('\nALL SANITY CHECKS PASSED ON ALL THREE MONTAGES.\n');
else
    fprintf('\n*** ONE OR MORE CHECKS FAILED -- do not use these baselines. ***\n');
end

function s = pf(x); if x, s='PASS'; else, s='**FAIL**'; end; end

function d = thesis_dir(varargin)
% Folder inside this repository
r = fileparts(mfilename('fullpath'));
while ~exist(fullfile(r,'README.md'),'file')
    p = fileparts(r);
    if isempty(p) || strcmp(p,r), error('thesis_dir:noRoot','no repo root'); end
    r = p;
end
d = [fullfile(r, varargin{:}) filesep];
end
