%A1_CONTACT_IMPEDANCE_SWEEP  Test how much contact impedance changes the readings and the difference signal.

% Aim served: Aim 1 (validation).
%
% Contact impedance is not known for this model, and two sources give values
% 14 times apart (Cheng 1989, McDermott 2020). This sweeps the range and checks
% how much it changes the absolute voltages and the difference signal dV/V.
%
% Inputs : a1_build_head_model
% Outputs: results/aim1/a1_contact_impedance_sweep.csv and .mat
% Imogen Ryan - Honours thesis.

clear; close all force;
addpath(fileparts(mfilename('fullpath')));
eidors_msg('log_level',1);
outd = thesis_dir('results','aim1');  if ~exist(outd,'dir'), mkdir(outd); end

Z_CHENG    = 5.5e-3;    % contact impedance (Ohm.m^2), Cheng 1989
Z_MCDERM   = 7.9e-2;    % contact impedance (Ohm.m^2), McDermott 2020
z_list = [1e-4 5e-4 1e-3 Z_CHENG 1e-2 3e-2 Z_MCDERM 1.5e-1 5e-1 1];

% Fixed perturbation: 10 mm diameter sphere, deep in the brain
SPH_CTR_mm = [0 0 20];  SPH_D_mm = 10;

fprintf('\n========= A1: CONTACT IMPEDANCE SWEEP =========\n');
fprintf('Cheng 1989 derived value  : %.4g Ohm.m^2\n', Z_CHENG);
fprintf('McDermott 2020 derived    : %.4g Ohm.m^2   (%.0fx larger)\n\n', Z_MCDERM, Z_MCDERM/Z_CHENG);
fprintf('%-12s %-13s %-13s %-14s %-13s\n', ...
    'z (Ohm.m^2)','RMS v (mV)','median|Z| Ohm','dV/V','vs Cheng dV/V');

T = table(); ref_dvv = NaN;
for z = z_list
    [img,info] = a1_build_head_model(struct('z_contact',z));
    v0 = fwd_solve(img); vb = v0.meas;

    % Same sphere at every contact impedance, so contact impedance is the only
    % thing changing. SPH_D_mm is a diameter, so it is halved for the radius.
    xyzr = [SPH_CTR_mm(:)/1000; (SPH_D_mm/2)/1000];
    map = mk_c2f_circ_mapping(img.fwd_model, xyzr);   % fraction of each element inside
    f = full(map(:,1)); f(~isfinite(f))=0; f = max(min(f,1),0);
    keep = false(size(f)); keep(info.brain_idx)=true; f(~keep)=0;   % brain only
    img2 = img; img2.elem_data = img.elem_data + f.*(info.sigma.csf - info.sigma.brain);

    % Relative size of the difference signal
    dvv = norm(fwd_solve(img2).meas - vb)/norm(vb);

    % Keep the result at the Cheng value, so every row can be given as a
    % percentage change from it
    if abs(z-Z_CHENG) < 1e-12, ref_dvv = dvv; end
    T = [T; table(z, rms(vb)*1e3, median(abs(vb))/info.opt.amp, dvv, ...
        'VariableNames',{'z_contact','rms_v_mV','medZ_Ohm','dVoverV'})];
end
T.dVV_rel_pct = 100*(T.dVoverV - ref_dvv)/ref_dvv;
for k = 1:height(T)
    fprintf('%-12.4g %-13.4f %-13.3f %-14.5e %+-13.2f\n', ...
        T.z_contact(k), T.rms_v_mV(k), T.medZ_Ohm(k), T.dVoverV(k), T.dVV_rel_pct(k));
end

iC = find(abs(T.z_contact-Z_CHENG)<1e-12);
iM = find(abs(T.z_contact-Z_MCDERM)<1e-12);
fprintf('\n--- ACROSS THE TWO SOURCED VALUES (Cheng -> McDermott, a 14x change) ---\n');
fprintf('  absolute RMS voltage changes by : %+.1f%%\n', ...
    100*(T.rms_v_mV(iM)-T.rms_v_mV(iC))/T.rms_v_mV(iC));
fprintf('  the DIFFERENCE signal dV/V by   : %+.2f%%\n', T.dVV_rel_pct(iM));
fprintf('\n  Across the FULL swept range (%.0e to %.0e, a %.0fx span):\n', ...
    min(z_list), max(z_list), max(z_list)/min(z_list));
fprintf('    absolute RMS voltage spans %.2fx\n', max(T.rms_v_mV)/min(T.rms_v_mV));
fprintf('    dV/V spans                 %.3fx\n', max(T.dVoverV)/min(T.dVoverV));

writetable(T,[outd 'a1_contact_impedance_sweep.csv']);
save([outd 'a1_contact_impedance_sweep.mat'],'T','Z_CHENG','Z_MCDERM');

fprintf('\nSAVED %s\n', [outd 'a1_contact_impedance_sweep.csv']);

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
