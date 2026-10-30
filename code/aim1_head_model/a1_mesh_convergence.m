%A1_MESH_CONVERGENCE  Check that the answer stops changing as the mesh gets finer.

% Aim served: Aim 1 (mesh convergence).
%
% Solves the healthy baseline and the change from a fixed 10 mm test sphere on
% all seven Tizzard adult meshes, which are one head at different densities.
% This shows numerical convergence only, not that one head represents everyone.
% The test sphere is a numerical check, not an EPVS simulation.
%
% Inputs : a1_build_head_model, all seven SAH*.mat meshes
% Outputs: results/aim1/a1_mesh_convergence.csv and .mat
% Imogen Ryan - Honours thesis.

clear; close all force;
addpath(fileparts(mfilename('fullpath')));
eidors_msg('log_level',1);   % quieten c2f progress messages
outd = thesis_dir('results','aim1');
if ~exist(outd,'dir'), mkdir(outd); end

meshes = {'SAH028','SAH031','SAH035','SAH049','SAH077','SAH132','SAH262'};
SPH_CTR_mm = [0 0 20];      % a fixed point inside the brain, in mm
SPH_D_mm   = 10;            % test sphere diameter (mm), resolved on all meshes
ELEC_AREA_m2 = 250e-6;      % target electrode area (m^2), the same on every mesh

% Pin every mesh to the same electrode centres, taken from the finest mesh,
% so electrode position does not change between meshes
[~, iref] = a1_build_head_model(struct('mesh',meshes{end}, ...
                'elec_target_area_m2',ELEC_AREA_m2));
ELEC_TARGETS = iref.elec_centre_m;

T = table();

fprintf('\n============== A1: MESH CONVERGENCE ==============\n');
fprintf('Electrodes pinned to SAH262 centres; electrode area fixed at %g mm^2.\n', ...
    ELEC_AREA_m2*1e6);
fprintf('Test sphere: %g mm diameter at [%g %g %g] mm\n\n', SPH_D_mm, SPH_CTR_mm);
fprintf('%-8s %-9s %-11s %-10s %-12s %-11s %-13s\n', ...
    'mesh','n_elem','elec mm^2','snap mm','RMS v (mV)','Z_1 (Ohm)','dV/V');

for m = 1:numel(meshes)
    t0 = tic;
    % Fixed electrode area across meshes, so that only mesh density changes
    [img,info] = a1_build_head_model(struct('mesh',meshes{m}, ...
                    'elec_target_area_m2', ELEC_AREA_m2, ...
                    'elec_targets', ELEC_TARGETS));
    vb = fwd_solve(img); v0 = vb.meas;

    % Insert the test sphere, [x; y; z; radius] in metres
    xyzr = [SPH_CTR_mm(:)/1000; (SPH_D_mm/2)/1000];

    % Fraction of each element inside the sphere, renormalised so the mapped
    % volume equals 4/3*pi*r^3. Fractions keep the same physical sphere on every
    % mesh, instead of rounding it to whole elements.
    map  = mk_c2f_circ_mapping(img.fwd_model, xyzr);
    f    = full(map(:,1));                 % fraction of each element inside the sphere
    f(~isfinite(f)) = 0; f = max(min(f,1),0);   % no NaN, limit to [0,1]

    % Restrict the perturbation to brain tissue
    inbrain = false(size(f)); inbrain(info.brain_idx) = true;
    f(~inbrain) = 0;

    % Volume perturbed, to show the same sphere was inserted on every mesh
    Vpert = sum(f .* info.elem_volume_m3);                % m^3

    % Linear volume mixing: an element half inside the sphere gets half the
    % CSF-brain difference. Used only because the test needs the same rule on every mesh.
    img2 = img;
    img2.elem_data = img.elem_data + f.*(info.sigma.csf - info.sigma.brain);
    dv = fwd_solve(img2).meas - v0;        % perturbed minus healthy

    % Z1 is the transfer impedance of the first channel (Ohm), a spot check on
    % the baseline. dvv is ||dv||/||v0||, the size of the whole difference
    % relative to the baseline, which is the quantity whose convergence matters.
    Z1 = abs(v0(1))/info.opt.amp;
    dvv = norm(dv)/norm(v0);
    tt = toc(t0);
    fprintf('%-8s %-9d %-11.1f %-10.2f %-12.4f %-11.4f %-13.5e\n', ...
        meshes{m}, size(img.fwd_model.elems,1), median(info.elec_area_m2)*1e6, ...
        mean(info.elec_snap_err_mm), rms(v0)*1e3, Z1, dvv);

    T = [T; table(string(meshes{m}), size(img.fwd_model.elems,1), ...
        median(info.elec_area_m2)*1e6, mean(info.elec_snap_err_mm), ...
        rms(v0)*1e3, Z1, dvv, Vpert*1e9, tt, ...
        'VariableNames',{'mesh','n_elems','elec_area_mm2','elec_snap_mm', ...
        'rms_v_mV','Z_first_Ohm','dVoverV','V_perturbed_mm3','solve_s'})];
end

% Convergence relative to the finest mesh
% There is no exact answer, so each mesh is compared with SAH262 as a
% percentage deviation. Small, shrinking deviations mean small discretisation error.
ref = height(T);
T.err_rms_pct  = 100*abs(T.rms_v_mV   - T.rms_v_mV(ref)) ./T.rms_v_mV(ref);
T.err_dVV_pct  = 100*abs(T.dVoverV    - T.dVoverV(ref))  ./T.dVoverV(ref);

fprintf('\n--- deviation from the finest mesh (SAH262) ---\n');
fprintf('%-8s %-9s %-16s %-16s\n','mesh','n_elem','baseline err %','dV/V err %');
for m=1:height(T)
    fprintf('%-8s %-9d %-16.2f %-16.2f\n', T.mesh(m), T.n_elems(m), ...
        T.err_rms_pct(m), T.err_dVV_pct(m));
end
fprintf('\nNOTE: the BASELINE converging does not mean the DIFFERENCE has converged.\n');
fprintf('The difference is the quantity the thesis measures, so its column is the\n');
fprintf('one that matters. No acceptance tolerance is proposed here -- a defensible\n');
fprintf('one should be tied to the Aim 3 detection threshold and agreed.\n');

writetable(T,[outd 'a1_mesh_convergence.csv']);
save([outd 'a1_mesh_convergence.mat'],'T','SPH_CTR_mm','SPH_D_mm');

fprintf('\nSAVED %s\n', [outd 'a1_mesh_convergence.csv']);

function d = thesis_dir(varargin)
% Folder inside this repository
r = fileparts(mfilename('fullpath'));
while ~exist(fullfile(r,'README.md'),'file')
    p = fileparts(r);
    if isempty(p) || strcmp(p,r)
        error('thesis_dir:noRoot', ...
            'Could not find the repository root (README.md) above %s', mfilename('fullpath'));
    end
    r = p;
end
d = [fullfile(r, varargin{:}) filesep];
end
