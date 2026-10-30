function sigma = a2_mix_conductivity(sigma_b, sigma_i, f, rule)

%A2_MIX_CONDUCTIVITY  Turn a fluid volume fraction into an effective conductivity.

% Aim served: Aim 2.
%
%   sigma = a2_mix_conductivity(sigma_b, sigma_i, f, rule)
%     sigma_b  background (brain) conductivity, S/m
%     sigma_i  inclusion (CSF) conductivity, S/m
%     f        fluid volume fraction per element, 0 to 1
%     rule     'linear' (L = 0, default) or 'maxwell_garnett' (L = 1/3, the main case)
% Imogen Ryan - Honours thesis.

if nargin < 4, rule = 'linear'; end

% A fraction outside [0,1] means an error in the caller, so stop rather than clamp.
if any(f < -1e-12) || any(f > 1+1e-12)
    error('a2_mix_conductivity:badFraction', ...
        'Volume fraction outside [0,1]: min %.3g, max %.3g', min(f), max(f));
end
f = max(min(f,1),0);            % rounding errors only

switch lower(rule)

    case 'linear'
        sigma = f.*sigma_i + (1-f).*sigma_b;

    case {'maxwell_garnett','mg'}
        beta  = (sigma_i - sigma_b) ./ (sigma_i + 2*sigma_b);
        denom = 1 - f.*beta;
        % denom reaches zero only if sigma_b <= 0. Checked so that an Inf
        % cannot pass into the solve.
        if any(denom <= 0)
            error('a2_mix_conductivity:singular', ...
                ['Maxwell-Garnett denominator non-positive. Check that ' ...
                 'sigma_b = %.4g is positive.'], sigma_b);
        end
        sigma = sigma_b .* (1 + 2*f.*beta) ./ denom;   % Maxwell Garnett 1904

    otherwise
        error('a2_mix_conductivity:rule', ...
            'Unknown mixing rule "%s". Use ''linear'' or ''maxwell_garnett''.', rule);
end
end
