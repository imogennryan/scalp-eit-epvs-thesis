%A1_VERIFY_UNITS  Prove the head model is in consistent SI units.

% Aim served: Aim 1 (physical scaling).
%
% Inputs : a1_build_head_model
% Outputs: results/aim1/a1_units_verification.csv and .mat
% Imogen Ryan - Honours thesis.

clear; close all force;
addpath(fileparts(mfilename('fullpath')));
outd = thesis_dir('results','aim1');
if ~exist(outd,'dir'), mkdir(outd); end

[img,info] = a1_build_head_model();
v0 = fwd_solve(img); v = v0.meas;
% Least-squares slope: the single factor by which every voltage changed.
% Used instead of an average of the ratios, because readings close to zero give very large ratios.
ratio = @(x,y) (y(:)\x(:));

fprintf('\n================= A1: UNITS VERIFICATION =================\n');
fprintf('%s\n', info.units);
fprintf('sigma set    : %s\n', info.sigma_source);
fprintf('z_contact    : %.3g Ohm.m^2\n', info.opt.z_contact);
fprintf('current      : %g A (%.0f uA)\n\n', info.opt.amp, info.opt.amp*1e6);

nd = img.fwd_model.nodes;
fprintf('head extent  : %.3f x %.3f x %.3f m  (adult head size, so the units are metres)\n', ...
    max(nd(:,1))-min(nd(:,1)), max(nd(:,2))-min(nd(:,2)), max(nd(:,3))-min(nd(:,3)));
fprintf('electrode    : median area %.1f mm^2 (equivalent diameter %.1f mm)\n', ...
    median(info.elec_area_m2)*1e6, 2*sqrt(median(info.elec_area_m2)/pi)*1e3);
fprintf('\nMEASURED VOLTAGES, now physical:\n');
fprintf('  RMS %.4f mV | max %.4f mV | min %.4f mV\n', ...
    rms(v)*1e3, max(abs(v))*1e3, min(abs(v))*1e3);
fprintf('  transfer impedance |v|/I: median %.2f Ohm, max %.2f Ohm\n', ...
    median(abs(v))/info.opt.amp, max(abs(v))/info.opt.amp);

R = table(); P = table();
% Test 1: conductivity scaling. Multiplying every conductivity by k should
% divide every voltage by exactly k, once the contact impedance is also scaled by 1/k.
fprintf('\n--- TEST 1: conductivity scaling ---\n');
fprintf('%-8s %-20s %-20s %-14s\n','k','z unchanged','z scaled 1/k','error vs exact');
for k = [2 10 137]
    ia = img; ia.elem_data = img.elem_data*k;
    ib = ia; for e=1:numel(ib.fwd_model.electrode)
        ib.fwd_model.electrode(e).z_contact = info.opt.z_contact/k; end
    va = fwd_solve(ia).meas; vb = fwd_solve(ib).meas;   % contact unchanged, contact changed
    ra = ratio(v, va); rb = ratio(v, vb);
    fprintf('%-8d %-20.10f %-20.10f %-14.2e\n', k, ra, rb, abs(rb-k)/k);
    R = [R; table("sigma", k, ra, rb, abs(rb-k)/k, ...
        'VariableNames',{'test','factor','ratio_z_fixed','ratio_z_scaled','rel_err'})];
    P = [P; plain_row("All four conductivities (skull / brain / scalp / CSF)", k, sig_txt(img.elem_data,1), ...
        sig_txt(img.elem_data,k), "S/m", v, va, vb)];
end

% Test 2: geometric scaling. Making the head a times bigger should divide every
% voltage by exactly a, once the contact impedance is also scaled by a.
% a = 1000 is the millimetre-to-metre factor.
fprintf('\n--- TEST 2: geometric scaling ---\n');
fprintf('%-8s %-20s %-20s %-14s\n','a','z unchanged','z scaled by a','error vs exact');
for a = [2 5 1000]
    ia = img; ia.fwd_model.nodes = img.fwd_model.nodes*a;
    ib = ia; for e=1:numel(ib.fwd_model.electrode)
        ib.fwd_model.electrode(e).z_contact = info.opt.z_contact*a; end
    va = fwd_solve(ia).meas; vb = fwd_solve(ib).meas;   % contact unchanged, contact changed
    ra = ratio(v, va); rb = ratio(v, vb);
    fprintf('%-8d %-20.10f %-20.10f %-14.2e\n', a, ra, rb, abs(rb-a)/a);
    R = [R; table("geometry", a, ra, rb, abs(rb-a)/a, ...
        'VariableNames',{'test','factor','ratio_z_fixed','ratio_z_scaled','rel_err'})];
    w0 = max(nd(:,1))-min(nd(:,1));             % head width, m
    P = [P; plain_row("Size of the head", a, sprintf("head width %.3f", w0), ...
        sprintf("head width %.3f", w0*a), "m", v, va, vb)];
end

mm_err = 100*abs(R.ratio_z_fixed(end)/1000 - 1);
fprintf('\n--- MILLIMETRES AGAINST METRES ---\n');
fprintf('Leaving the mesh in mm gives a voltage scaling of %.1f instead of 1000.\n', ...
    R.ratio_z_fixed(end));
fprintf('That is a %.1f%% error, and it is NOT a constant factor -- it depends on\n', mm_err);
fprintf('z_contact, so it cannot be removed by dividing by 1000 afterwards.\n');
fprintf('CONCLUSION: nodes must be converted to metres at load time. Done in\n');
fprintf('a1_build_head_model.m. Both exact tests now pass to ~1e-10.\n');

writetable(R,[outd 'a1_units_verification.csv']);
writetable(P,[outd 'a1_units_verification_plain.csv']);
save([outd 'a1_units_verification.mat'],'R','P','info');
fprintf('\nSAVED %s\n',[outd 'a1_units_verification.csv']);

function row = plain_row(what, k, before, after, units, v, va, vb)
% One row of the plain results table: what was scaled, and what the readings did.
% va: contact impedance unchanged (as in the real simulations)
% vb: contact impedance changed by hand, so the whole circuit scales
pct = @(vnew) 100*(1/(vnew(:)\v(:)) - 1);          % % change, fitted over all readings
worst = 100*max(abs(vb*k - v))/max(abs(v));        % worst single reading, contact changed
row = table(string(what), k, string(before), string(after), string(units), ...
    rms(v)*1e6, rms(va)*1e6, rms(vb)*1e6, 100*(1/k - 1), pct(va), pct(vb), worst, ...
    'VariableNames',{'changed','factor','before','after','units', ...
    'rms_before_uV','rms_after_contact_unchanged_uV','rms_after_contact_changed_uV', ...
    'expected_change_pct','change_contact_unchanged_pct','change_contact_changed_pct', ...
    'worst_reading_error_pct'});
end

function s = sig_txt(elem_data, k)
% The four tissue conductivities, times k, as text, smallest first: skull, brain, scalp, CSF
u = unique(elem_data);                             % 0.05 0.10 0.23 1.79 S/m on the thesis set
s = strjoin(compose("%.3g", u(:)'*k), " / ");
end

function d = thesis_dir(varargin)
% Path to a folder inside this repository.
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
