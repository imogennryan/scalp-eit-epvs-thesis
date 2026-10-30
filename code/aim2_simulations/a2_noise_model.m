function [noise, ninfo] = a2_noise_model(v_baseline, snr_db, n_real, seed)

%A2_NOISE_MODEL  Gaussian instrument noise scaled to the healthy baseline frame.

% Aim served: Aim 2.
%
%   [noise, ninfo] = a2_noise_model(v_baseline, snr_db, n_real, seed)
%     v_baseline  healthy baseline measurement vector, V
%     snr_db      instrument SNR, dB (scalar or vector)
%     n_real      number of noise realisations (default 100)
%     seed        random seed (default 0)
%     noise       numel(v_baseline) x n_real x numel(snr_db), V
%     ninfo       struct recording the noise settings used
%
% The noise is scaled to the baseline frame, not to the difference, because an
% instrument's noise floor does not depend on what is being measured
% (McDermott 2020). SNR is an amplitude ratio, 10^(dB/20) (Mamatjan 2011).
% Noise is the same on every channel and frame averaging is not modelled.
% Imogen Ryan - Honours thesis.

if nargin < 3 || isempty(n_real), n_real = 100; end
if nargin < 4 || isempty(seed),   seed   = 0;   end

v_baseline = v_baseline(:);
nm  = numel(v_baseline);
snr_db = snr_db(:)';
ns  = numel(snr_db);

% Fixed seed so the noise is reproducible.
rs = RandStream('mt19937ar','Seed',seed);

% dB to linear amplitude ratio.
snr_lin = 10.^(snr_db/20);

% Noise is referenced to the baseline frame, not the difference.
ref_norm = norm(v_baseline);

noise = zeros(nm, n_real, ns);
for s = 1:ns
    for k = 1:n_real
        w = randn(rs, nm, 1);
        % Scale so that norm(v_baseline)/norm(w) equals snr_lin, as in EIDORS add_noise.
        noise(:,k,s) = w * ref_norm / norm(w) / snr_lin(s);
    end
end

% Noise settings, saved with any result that uses them.
ninfo = struct();
ninfo.snr_db          = snr_db;
ninfo.snr_linear      = snr_lin;
ninfo.db_convention   = 'amplitude: SNR_linear = 10^(dB/20)';
ninfo.reference       = 'baseline measurement frame, not the difference';
ninfo.reference_norm  = ref_norm;
ninfo.n_realisations  = n_real;
ninfo.seed            = seed;
ninfo.channel_model   = 'uniform across channels; direction of bias unknown, see header';
ninfo.averaging       = 'not modelled; real systems gain sqrt(N) by frame averaging';
ninfo.units           = 'volts';
end
