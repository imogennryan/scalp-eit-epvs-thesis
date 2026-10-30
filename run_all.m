%RUN_ALL  Run every thesis script in order, from a clean workspace.

% Aim served: Aim 1 and Aim 2.
%
% Inputs : EIDORS 3.12 on the MATLAB path, models/tizzard_meshes/
% Outputs: models/, data/, results/ and figures/ inside this folder
% Imogen Ryan - Honours thesis.

assert(exist('eidors_obj','file')==2, ...
    'EIDORS is not on the MATLAB path. Run EIDORS startup.m first.');

here = fileparts(mfilename('fullpath'));

% Aim 1: build and validate the head model
aim1 = {'a1_build_montages'               % the three electrode layouts, needed by all later steps
        'a1_verify_mesh_family'           % Test 1, mesh quality
        'a1_geometry_verification'        % Tests 2 and 3, nesting and head size
        'a1_verify_units'                 % Test 4, scaling
        'a1_physical_plausibility'        % Test 5, power and current by tissue
        'a1_reciprocity_table'            % Test 6, reciprocity
        'a1_mesh_convergence'             % Test 7, mesh convergence
        'a1_contact_impedance_sweep'      % contact impedance test
        'a1_analytic_sphere_validation'   % Test 8, volume law and levelling off
        'a1_generate_baseline'            % healthy head for each layout
        'a1_methods_figures'              % Figures 3.1 to 3.3
        'a1_volume_law_figure'            % Figure 4.1
        'a1_saturation_figure'};          % Figure 4.2

% Aim 2: simulate the readings
aim2 = {'a2_run_sweep'                    % all fluid simulations
        'a2_add_noise'                    % noise and detection threshold
        'a2_size_sweep'                   % region size test
        'a2_results_summary'              % Results tables
        'a2_methods_figures'              % Figures 3.4 to 3.6
        'a2_methods_figure_csf'           % Figure 3.7
        'a2_noise_method_figure'          % Figure 3.8
        'a2_fluid_depth_position_figure'  % Figures 4.3 and 4.4
        'a2_margin_figure'                % Figure 4.5
        'a2_results_figures_v3'           % Figure 4.6
        'a2_sensitivity_figure'};         % Figure 4.7

steps = [fullfile(here,'code','aim1_head_model',strcat(aim1,'.m'));
         fullfile(here,'code','aim2_simulations',strcat(aim2,'.m'))];

t0 = tic;
for k = 1:numel(steps)
    [~,name] = fileparts(steps{k});
    fprintf('\n>>> [%d/%d] %s\n', k, numel(steps), name);
    tk = tic;
    run_one(steps{k});
    fprintf('    (%.1f s)\n', toc(tk));
end
fprintf('\nAll %d scripts finished in %.1f min.\n', numel(steps), toc(t0)/60);

function run_one(p)
% Each script runs in its own workspace, so its clear does not stop this loop
run(p);
end
