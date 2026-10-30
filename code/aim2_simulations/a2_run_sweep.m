%A2_RUN_SWEEP  Generate the full simulated dataset of voltage differences from the healthy baseline.

% Aim served: Aim 2.
%
% Axes: depth (25 to 60 mm on a radial line from a scalp point), scalp point,
% fluid volume, scenario, electrode montage and mixing rule
% ('maxwell_garnett', L = 1/3, main case; 'linear', L = 0, best-case bound).
% Scenario B (main case) is a 20 mm region holding a fluid volume fraction.
% Scenario A is one CSF-filled space of 3 or 5 mm. It is smaller than the
% mesh elements, so it is biased toward detection and is a bound only.
% Combinations deeper than the head allows from a scalp point are recorded
% as unreachable and not run. Noise is added later, in a2_add_noise.
% Each configuration has a key built from its parameters, so a stopped run
% resumes where it left off.
%
% Inputs : a1_build_head_model, a2_place_perturbation, a2_mix_conductivity,
%          a2_reachable_depth
% Outputs: data/aim2/a2_sweep_<montage>.mat      difference vectors + baseline
%          results/aim2/a2_sweep_index.csv       one row per configuration
%          results/aim2/a2_sweep_reachability.csv  max depth per position
%          results/aim2/a2_sweep_unreachable.csv   configurations that cannot exist
% Imogen Ryan - Honours thesis.

clear; close all force;
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(fileparts(here),'aim1_head_model'));
eidors_msg('log_level',1);
rng('default');                                   % fixed seed

datd = thesis_dir('data','aim2');    if ~exist(datd,'dir'), mkdir(datd); end
outd = thesis_dir('results','aim2'); if ~exist(outd,'dir'), mkdir(outd); end

%% Axis definitions
DEPTHS   = 25:5:60;
ANCHORS  = struct( ...
    'name', {'Cz','Fz','C3','Pz','T7','O1'}, ...
    'note', {'vertex, centrum semiovale territory', ... % Banerjee et al. 2017; Kim et al. 2021
             'frontal white matter', ...                % Braun et al. 2024; Piantino et al. 2021
             'left central, lateral line toward deep grey', ...
             'parietal, posterior coverage', ...
             'left temporal, most lateral line', ...
             'left occipital, most posterior line'});
BURDENS  = [0.0025 0.005 0.01 0.02 0.05 0.10];
MONTAGES = {'m19_1020','m32_1010','m16_twoplane'};
RULES    = {'linear','maxwell_garnett'};
REGION_D_MM = 20;                                 % Scenario B region diameter (mm)
SINGLE_D_MM = [3 5];                              % Scenario A space diameters (mm), Wardlaw
SAVE_EVERY  = 100;

% Look up each anchor's position by name in the saved 10-10 table. Anchors are
% physical scalp points, so they are the same for every montage.
Lm = load([thesis_dir('models','montages') 'a1_montages.mat'],'pos','names');
for a = 1:numel(ANCHORS)
    j = find(strcmp(Lm.names, ANCHORS(a).name), 1);
    assert(~isempty(j),'Anchor "%s" is not in the 10-10 position table.', ANCHORS(a).name);
    ANCHORS(a).xyz_m = Lm.pos(j,:)/1000;
end

fprintf('\n############### A2: FULL SWEEP ###############\n');

%% Reachable depth from each scalp position
% The mesh is the same for every montage, so this is computed once.
fprintf('\n--- reachable depth by scalp position ---\n');
[img0, info0] = a1_build_head_model();
R = table();
fprintf('%-6s %-14s %-34s\n','anchor','max depth mm','deepest usable centre (region / single)');
for a = 1:numel(ANCHORS)
    dm = a2_reachable_depth(img0, info0, ANCHORS(a).xyz_m);
    ANCHORS(a).dmax = dm;
    fprintf('%-6s %-14.1f %.1f mm / %.1f mm\n', ANCHORS(a).name, dm, ...
        dm - REGION_D_MM/2, dm - max(SINGLE_D_MM)/2);
    R = [R; table(string(ANCHORS(a).name), dm, dm-REGION_D_MM/2, ...
        dm-max(SINGLE_D_MM)/2, string(ANCHORS(a).note), ...
        'VariableNames',{'anchor','max_depth_mm','max_centre_region_mm', ...
        'max_centre_single_mm','note'})]; %#ok<AGROW>
end
writetable(R, [outd 'a2_sweep_reachability.csv']);
fprintf('\nMax depth varies %.1f to %.1f mm by position. This is a property of\n', ...
    min([ANCHORS.dmax]), max([ANCHORS.dmax]));
fprintf('head shape, not of the model, and it limits the grid below.\n');

%% Build the configuration list
% A configuration is run only if the whole perturbation fits within the
% reachable depth. The others are recorded separately.
C = table(); U = table();
for mo = 1:numel(MONTAGES)
  for an = 1:numel(ANCHORS)
    for d = DEPTHS
      specs = [repmat({'B'},1,numel(BURDENS)), repmat({'A'},1,numel(SINGLE_D_MM))];
      sizes = [repmat(REGION_D_MM,1,numel(BURDENS)), SINGLE_D_MM];
      burds = [BURDENS, ones(1,numel(SINGLE_D_MM))];
      for q = 1:numel(specs)
        fits = ANCHORS(an).dmax >= d + sizes(q)/2;
        for r = 1:numel(RULES)
          rw = mkrow(specs{q}, MONTAGES{mo}, ANCHORS(an), d, sizes(q), burds(q), RULES{r});
          if fits, C = [C; rw]; else, U = [U; rw]; end %#ok<AGROW>
        end
      end
    end
  end
end
writetable(U, [outd 'a2_sweep_unreachable.csv']);

fprintf('\n%d configurations to run, %d unreachable (recorded, not skipped silently)\n', ...
    height(C), height(U));
fprintf('  scenario B (regional burden): %d\n', sum(C.scenario=="B"));
fprintf('  scenario A (single space)   : %d\n', sum(C.scenario=="A"));
% Depth range that every position can reach in scenario B.
core_depth = min([ANCHORS.dmax]) - REGION_D_MM/2;
fprintf('  balanced core for scenario B: depths up to %.0f mm at ALL %d positions\n\n', ...
    floor(core_depth/5)*5, numel(ANCHORS));

%% Resume from stored configurations
idxfile = [outd 'a2_sweep_index.csv'];
IDX = table(); done_keys = strings(0,1);
if exist(idxfile,'file')
    IDX = readtable(idxfile,'TextType','string');
    done_keys = IDX.key;
    fprintf('RESUMING: %d configurations already stored.\n', numel(done_keys));
end
todo = C(~ismember(C.key, done_keys), :);
if isempty(todo)
    fprintf('Dataset already complete. Delete %s to force a re-run.\n', idxfile); return
end
fprintf('%d configurations still to run.\n\n', height(todo));

%% Placements, computed once
% Placement depends only on anchor, depth and size, not on montage. Sharing it
% means every montage sees the identical perturbation.
fprintf('--- computing placements (shared across montages) ---\n');
[Gp, pa, pd, ps] = findgroups(todo.anchor, todo.depth_mm, todo.size_mm);
PL = containers.Map('KeyType','char','ValueType','any');
t0 = tic;
for g = 1:max(Gp)
    ai = find(strcmp({ANCHORS.name}, pa(g)), 1);
    o = struct('mode','radial','depth_mm',pd(g),'diam_mm',ps(g), ...
               'anchor_xyz_m',ANCHORS(ai).xyz_m);
    [frac, place] = a2_place_perturbation(img0, info0, o);
    % Stored sparse, as only a few hundred elements are non-zero.
    PL(plkey(pa(g),pd(g),ps(g))) = struct('frac',sparse(frac),'place',place);
end
fprintf('  %d placements in %.1f s\n\n', max(Gp), toc(t0));

%% Run
t_start = tic; n_done = 0;
for mo = 1:numel(MONTAGES)
    sub_mo = todo(todo.montage==string(MONTAGES{mo}), :);
    if isempty(sub_mo), continue; end

    % Baseline solved on the identical model for each montage.
    [img, info] = a1_build_head_model(struct('montage', MONTAGES{mo}));
    vb = fwd_solve(img).meas;
    fprintf('=== %s: %d electrodes, %d measurements, baseline RMS %.4f mV ===\n', ...
        MONTAGES{mo}, numel(img.fwd_model.electrode), numel(vb), rms(vb)*1e3);

    dvfile = [datd sprintf('a2_sweep_%s.mat', MONTAGES{mo})];
    if exist(dvfile,'file'), S = load(dvfile); DV = S.DV;
    else, DV = containers.Map('KeyType','char','ValueType','any'); end

    for k = 1:height(sub_mo)
        cfg = sub_mo(k,:);
        P   = PL(plkey(cfg.anchor, cfg.depth_mm, cfg.size_mm));
        frac = full(P.frac); place = P.place;

        fluid = frac * cfg.burden_frac;
        m = fluid > 0;
        img2 = img;
        img2.elem_data(m) = a2_mix_conductivity( ...
            info.sigma.brain, info.sigma.csf, fluid(m), char(cfg.rule));

        dv = fwd_solve(img2).meas - vb;
        DV(char(cfg.key)) = dv;                    % keyed by configuration

        [pk, ipk] = max(abs(dv));
        IDX = [IDX; [cfg, table( ...
            place.depth_mm_actual, place.vol_mm3, place.vol_fraction_delivered, ...
            place.fits_in_brain, place.n_elems_touched, ...
            cfg.burden_frac*place.vol_mm3, numel(vb), rms(vb)*1e3, ...
            norm(dv)/norm(vb), pk*1e3, pk/abs(vb(ipk)), ipk, ...
            place.centre_mm(1), place.centre_mm(2), place.centre_mm(3), ...
            'VariableNames',{'depth_actual_mm','pert_vol_mm3','vol_frac_delivered', ...
            'fits_in_brain','n_elems','fluid_vol_mm3','n_meas','baseline_rms_mV', ...
            'dVoverV','peak_dV_mV','peak_dV_rel','peak_channel', ...
            'cx_mm','cy_mm','cz_mm'})]]; %#ok<AGROW>

        n_done = n_done + 1;
        if mod(n_done, SAVE_EVERY) == 0
            save(dvfile,'DV','vb','-v7.3'); writetable(IDX, idxfile);
            fprintf('  ... %d/%d (%.0f s), saved\n', n_done, height(todo), toc(t_start));
        end
    end
    save(dvfile,'DV','vb','-v7.3'); writetable(IDX, idxfile);
    fprintf('  %s complete\n', MONTAGES{mo});
end
writetable(IDX, idxfile);

fprintf('\n############### SWEEP COMPLETE ###############\n');
fprintf('%d configurations in %.1f s (%.2f s each)\n', n_done, toc(t_start), toc(t_start)/max(n_done,1));

%% Dataset checks
fprintf('\n=============== DATASET INTEGRITY ===============\n');
c1 = height(IDX) == height(C);
fprintf('  every runnable configuration present : %s (%d of %d)\n', pf(c1), height(IDX), height(C));
c2 = numel(unique(IDX.key)) == height(IDX);
fprintf('  no duplicate configuration keys      : %s\n', pf(c2));
c3 = all(isfinite(IDX.dVoverV)) && all(IDX.dVoverV > 0);
fprintf('  all responses finite and non-zero    : %s\n', pf(c3));
c4 = all(IDX.fits_in_brain == 1);
fprintf('  every perturbation fully in brain    : %s (%d clipped)\n', pf(c4), sum(IDX.fits_in_brain==0));
% The placement is shared, so the volume must be identical across montages.
c5 = true;
for s = ["B","A"]
    g = IDX(IDX.scenario==s,:);
    [GG,~,~,~] = findgroups(g.anchor, g.depth_mm, g.size_mm);
    for q = 1:max(GG)
        v = g.pert_vol_mm3(GG==q);
        c5 = c5 && (max(v)-min(v)) < 1e-6;
    end
end
fprintf('  identical perturbation across montages: %s\n', pf(c5));
ok = c1 && c2 && c3 && c4 && c5;
fprintf('\n  ALL CHECKS PASS: %s\n', pf(ok));
if ~ok, fprintf('  Do not analyse this dataset until the failing check is understood.\n'); end
fprintf('\nindex -> %s\ndata  -> %s\n', idxfile, datd);

%% Functions
function r = mkrow(scen, montage, anchor, depth, size_mm, burden, rule)
% Key built from the parameters that define the configuration.
key = sprintf('%s|%s|%s|d%g|s%g|b%g|%s', scen, montage, anchor.name, ...
              depth, size_mm, burden, rule);
r = table(string(key), string(scen), string(montage), string(anchor.name), ...
    anchor.xyz_m(1), anchor.xyz_m(2), anchor.xyz_m(3), anchor.dmax, ...
    depth, size_mm, burden, string(rule), string(anchor.note), ...
    'VariableNames',{'key','scenario','montage','anchor', ...
    'anchor_x_m','anchor_y_m','anchor_z_m','anchor_dmax_mm', ...
    'depth_mm','size_mm','burden_frac','rule','anchor_note'});
end

function k = plkey(anchor, depth, size_mm)
k = sprintf('%s|d%g|s%g', char(anchor), depth, size_mm);
end

function s = pf(x); if x, s='PASS'; else, s='**FAIL**'; end; end

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
