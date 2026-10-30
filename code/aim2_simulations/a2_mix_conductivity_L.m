function sigma = a2_mix_conductivity_L(sigma_b, sigma_i, f, L)

%A2_MIX_CONDUCTIVITY_L  Effective conductivity as a function of inclusion shape.

% Aim served: Aim 2. Replaces an unsourced binary choice with one physical
%             parameter that can be argued about.
%
%   sigma = a2_mix_conductivity_L(sigma_b, sigma_i, f, L)
%     sigma_b  background (brain) conductivity, S/m
%     sigma_i  inclusion (CSF) conductivity, S/m
%     f        fluid volume fraction per element, 0 to 1
%     L        depolarisation factor of the inclusions, 0 to 1
% Imogen Ryan - Honours thesis.

if any(f < -1e-12) || any(f > 1+1e-12)
    error('a2_mix_conductivity_L:badFraction', ...
        'Volume fraction outside [0,1]: min %.3g, max %.3g', min(f), max(f));
end
if L < 0 || L > 1
    error('a2_mix_conductivity_L:badL', 'L must lie in [0,1], got %.3g', L);
end
f = max(min(f,1),0);

d = sigma_i - sigma_b;
denom = sigma_b + L .* (1-f) .* d;

% denom vanishes only for sigma_b <= 0 or L outside [0,1], both excluded above.
% Guarded anyway: a silent Inf would propagate into the solve and produce
% nonsense voltages instead of an error.

if any(abs(denom) < eps(sigma_b)*1e3)
    error('a2_mix_conductivity_L:singular', 'Degenerate denominator at L=%.3g.', L);
end

sigma = sigma_b .* (1 + f .* d ./ denom);
end
