%A1_BUILD_MONTAGES  Build, check and save the three EEG-standard electrode montages.

% Aim served: Aim 1 (electrode model) and Aim 2 (electrode configuration).
%
% The 19-electrode montage is the classical 10-20 set. The 32- and 16-electrode
% montages are built from the 10-10 set to match the electrode count and layout
% described by McDermott et al. 2020 and Ouypornkochagorn et al. 2023. Their
% exact positions are not published, so these are close matches, not copies.
% Electrode area snaps to whole boundary faces, so the 10 mm diameter target
% is only approximately reached.
%
% Inputs : a1_make_1010_positions, a1_build_head_model
% Outputs: models/montages/a1_montages.mat and results/aim1/a1_montages.csv
% Imogen Ryan - Honours thesis.

clear; close all force;
addpath(fileparts(mfilename('fullpath')));
eidors_msg('log_level',1);
outd = thesis_dir('results','aim1');      if ~exist(outd,'dir'), mkdir(outd); end
figd = thesis_dir('figures','aim1');      if ~exist(figd,'dir'), mkdir(figd); end   % for the figures drawn by later scripts
mdld = thesis_dir('models','montages');   if ~exist(mdld,'dir'), mkdir(mdld); end

MESH = 'SAH262';
[pos, names] = a1_make_1010_positions(MESH);

% Classical 10-20 set, 19 electrodes (international standard)
m19 = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8', ...
       'P7','P3','Pz','P4','P8','O1','O2'};

% 32-electrode 10-10 set, matched to McDermott 2020: 4 midline electrodes
% (Fz, Cz, Pz, FCz) and 14 symmetric pairs
m32 = [m19, {'AF7','AF8','FC5','FC1','FC2','FC6','CP5','CP1','CP2','CP6', ...
             'PO7','PO8','FCz'}];

% 16-electrode two-plane set, matched in count to Ouypornkochagorn 2023
% (a lower plane of 10, an upper plane of 6)
m16_lower = {'Fp1','Fp2','F7','F8','T7','T8','P7','P8','O1','O2'};   % 10
m16_upper = {'F3','F4','C3','C4','P3','P4'};                          % 6
m16 = [m16_lower, m16_upper];

M = struct('name',{'m19_1020','m32_1010','m16_twoplane'}, ...
           'labels',{m19,m32,m16}, ...
           'source',{'classical 10-20, international standard', ...
                     'count-matched to McDermott 2020, constructed from 10-10', ...
                     'count- and structure-matched to Ouypornkochagorn 2023'});

fprintf('\n============ A1: ELECTRODE MONTAGES (%s) ============\n', MESH);
T = table();
for k = 1:numel(M)
    lab = M(k).labels;
    % Stop on a duplicated or unknown electrode label
    assert(numel(unique(lab))==numel(lab), 'duplicate label in %s', M(k).name);
    [tf,loc] = ismember(lab, names);       % loc = row of each label in the 10-10 list
    assert(all(tf), 'unknown label(s) in %s: %s', M(k).name, strjoin(lab(~tf),', '));
    P = pos(loc,:);                        % this montage's positions, in mm
    M(k).pos_mm = P;

    % Build a model with this montage to find the electrode area that can be
    % reached, since electrodes snap to whole boundary faces, and how far each
    % electrode moved from its ideal 10-10 position.
    [img,info] = a1_build_head_model(struct('mesh',MESH, ...
        'elec_targets', P/1000, 'elec_target_area_m2', 78.5e-6));
    v = fwd_solve(img);
    M(k).elec_area_mm2 = info.elec_area_m2*1e6;
    M(k).snap_mm       = info.elec_snap_err_mm;
    nmeas = numel(v.meas);

    fprintf('\n%-14s %d electrodes | %s\n', M(k).name, numel(lab), M(k).source);
    fprintf('   measurements (adjacent, no_meas_current) : %d\n', nmeas);
    fprintf('   electrode area  : median %.1f mm^2 (target 78.5, i.e. 10 mm dia)\n', ...
        median(M(k).elec_area_mm2));
    fprintf('   equivalent diam : median %.1f mm\n', 2*sqrt(median(M(k).elec_area_mm2)/pi));
    fprintf('   snap to surface : median %.2f mm, max %.2f mm\n', ...
        median(M(k).snap_mm), max(M(k).snap_mm));
    fprintf('   baseline RMS    : %.4f mV\n', rms(v.meas)*1e3);
    T = [T; table(string(M(k).name), numel(lab), nmeas, median(M(k).elec_area_mm2), ...
        2*sqrt(median(M(k).elec_area_mm2)/pi), median(M(k).snap_mm), max(M(k).snap_mm), ...
        rms(v.meas)*1e3, string(M(k).source), 'VariableNames', ...
        {'montage','n_elec','n_meas','elec_area_mm2','elec_diam_mm', ...
         'snap_med_mm','snap_max_mm','rms_v_mV','provenance'})];
end

fprintf('\nNOTE: electrode area is quantised to whole boundary faces, so the 10 mm\n');
fprintf('target from Ouypornkochagorn 2023 is not exactly achievable on this mesh.\n');

save([mdld 'a1_montages.mat'],'M','pos','names','MESH');
writetable(T,[outd 'a1_montages.csv']);

fprintf('\nSAVED %s\nSAVED %s\n', [mdld 'a1_montages.mat'], [outd 'a1_montages.csv']);

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
