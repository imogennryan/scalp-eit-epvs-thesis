%A2_SIZE_SWEEP  Test whether accumulation size matters once fluid volume is held fixed.

% Aim served: Aim 2.
%
% The same fluid volume is spread through regions of different diameter, at
% three values of the shape factor L. Everything else is held fixed. The
% burden is set from the delivered region volume, so the fluid volume is
% matched even where the region is clipped at the brain boundary.
%
% Inputs : a1_build_head_model, a2_place_perturbation, a2_mix_conductivity_L,
%          a2_mix_conductivity
% Outputs: results/aim2/a2_size_sweep.csv
% Imogen Ryan - Honours thesis.

clear; close all force;
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(fileparts(here),'aim1_head_model'));
eidors_msg('log_level',1);
rng('default');

outd = thesis_dir('results','aim2'); if ~exist(outd,'dir'), mkdir(outd); end

DIAMS   = [12 16 20 26 32];        % region diameter (mm)
V_FLUID = [100 200 400];           % fluid volume held constant across DIAMS (mm^3)
DEPTHS  = [30 40];                 % depth below the scalp point (mm)
ANCHORS = {'Cz','Fz'};             % scalp points above Alzheimer's and brain injury regions
LVALS   = [0 1/3 0.5];             % shape factor L: tubes along the current, spheres, tubes across

fprintf('\n############### A2: MATCHED-VOLUME SIZE SWEEP ###############\n');

%% Check: the general formula reproduces both earlier mixing rules
% If it does not, the earlier results and this one are inconsistent, so the run stops.
ftest = (0:0.1:1)';
e0 = max(abs(a2_mix_conductivity_L(0.10,1.79,ftest,0)   - a2_mix_conductivity(0.10,1.79,ftest,'linear')));
e3 = max(abs(a2_mix_conductivity_L(0.10,1.79,ftest,1/3) - a2_mix_conductivity(0.10,1.79,ftest,'maxwell_garnett')));
fprintf('\n--- equivalence with the existing rules ---\n');
fprintf('  L=0   vs linear          : max diff %.2e S/m\n', e0);
fprintf('  L=1/3 vs maxwell_garnett : max diff %.2e S/m\n', e3);
assert(e0 < 1e-12 && e3 < 1e-12, 'Generalised mixing does not reproduce the existing rules.');
fprintf('  PASS -- the two old rules are L=0 and L=1/3 of one formula.\n');

%% Model and baseline
[img, info] = a1_build_head_model();
vb = fwd_solve(img).meas;
fprintf('\nbaseline: %d measurements, RMS %.4f mV\n', numel(vb), rms(vb)*1e3);

Lm = load([thesis_dir('models','montages') 'a1_montages.mat'],'pos','names');
T = table(); nrun = 0; t0 = tic;

for ia = 1:numel(ANCHORS)
  j = find(strcmp(Lm.names, ANCHORS{ia}),1);
  anchor = Lm.pos(j,:)/1000;
  for dp = DEPTHS
    for D = DIAMS
      o = struct('mode','radial','depth_mm',dp,'diam_mm',D,'anchor_xyz_m',anchor);
      [frac, place] = a2_place_perturbation(img, info, o);

      % Skip regions that do not fit in the brain, rather than reporting a
      % clipped region as a full one
      if place.vol_fraction_delivered < 0.95
        fprintf('  skip %s d%d D%d : only %.1f%% of the region fits\n', ...
            ANCHORS{ia}, dp, D, 100*place.vol_fraction_delivered);
        continue
      end

      for Vf = V_FLUID
        b = Vf / place.vol_mm3;                 % burden set to match the fluid volume
        if b > 1, continue, end                 % too much fluid for the region
        fluid = frac * b;  m = fluid > 0;
        vol_actual = sum(fluid .* info.elem_volume_m3) * 1e9;

        for L = LVALS
          img2 = img;
          img2.elem_data(m) = a2_mix_conductivity_L(info.sigma.brain, info.sigma.csf, fluid(m), L);
          dv = fwd_solve(img2).meas - vb;
          nrun = nrun + 1;

          T = [T; table(string(ANCHORS{ia}), dp, place.depth_mm_actual, D, ...
               place.vol_mm3, Vf, vol_actual, b, L, place.n_elems_touched, ...
               max(abs(dv))*1e3, norm(dv), rms(vb)*1e3, ...
               'VariableNames', {'anchor','depth_target_mm','depth_actual_mm', ...
               'diam_mm','region_vol_mm3','fluid_vol_target_mm3','fluid_vol_actual_mm3', ...
               'burden_frac','L','n_elems','peak_dV_mV','norm_dV','baseline_rms_mV'})]; %#ok<AGROW>
        end
      end
    end
  end
end
fprintf('\n%d solves in %.1f s\n', nrun, toc(t0));

%% Check the fluid volume match held
relerr = abs(T.fluid_vol_actual_mm3 - T.fluid_vol_target_mm3) ./ T.fluid_vol_target_mm3;
fprintf('\n--- did the matched-volume construction work? ---\n');
fprintf('  worst relative error in delivered fluid volume: %.2e\n', max(relerr));
assert(max(relerr) < 1e-6, 'Fluid volume was not matched; the comparison is invalid.');
fprintf('  PASS -- fluid volume is matched across diameters to better than 1e-6.\n');

writetable(T, [outd 'a2_size_sweep.csv']);

%% Results table
fprintf('\n=========== DOES SIZE MATTER AT FIXED FLUID VOLUME? ===========\n');
fprintf('peak dV in mV. Each row holds fluid volume constant and spreads it\n');
fprintf('through a region of increasing diameter.\n\n');
for L = LVALS
  fprintf('L = %.3f\n', L);
  fprintf('  %-8s %-8s', 'anchor', 'depth');
  fprintf('%-10s', 'Vf(mm3)');
  for D = DIAMS, fprintf('D=%-8d', D); end
  fprintf('%-10s\n','spread');
  for ia = 1:numel(ANCHORS)
    for dp = DEPTHS
      for Vf = V_FLUID
        k = T.anchor==string(ANCHORS{ia}) & T.depth_target_mm==dp & ...
            T.fluid_vol_target_mm3==Vf & abs(T.L-L)<1e-9;
        if sum(k) < 2, continue, end
        pk = T.peak_dV_mV(k); dd = T.diam_mm(k);
        [dd,ix] = sort(dd); pk = pk(ix);
        fprintf('  %-8s %-8d%-10d', ANCHORS{ia}, dp, Vf);
        vals = nan(1,numel(DIAMS));
        for q = 1:numel(dd), vals(DIAMS==dd(q)) = pk(q); end
        for q = 1:numel(DIAMS)
          if isnan(vals(q)), fprintf('%-10s','-'); else, fprintf('%-10.5f',vals(q)); end
        end
        fprintf('%+.1f dB\n', 20*log10(max(pk)/min(pk)));
      end
    end
  end
  fprintf('\n');
end

fprintf('SAVED %s   (%d rows)\n', [outd 'a2_size_sweep.csv'], height(T));

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
