function [dmax, prof] = a2_reachable_depth(img, info, anchor_xyz_m)

%A2_REACHABLE_DEPTH  Greatest depth reachable along a radial line inward from a scalp point.

% Aim served: Aim 2.
%
%   [dmax, prof] = a2_reachable_depth(img, info, anchor_xyz_m)
%     img           head model image
%     info          model information (not used)
%     anchor_xyz_m  scalp point the line starts from, m
%     dmax          greatest depth below the scalp along the line, mm
%     prof          struct with the depth profile along the line
%
% Depth is the distance to the nearest scalp point. It rises and then falls
% along the line, so the maximum depends on where on the scalp the line starts.
% Imogen Ryan - Honours thesis.

nodes  = img.fwd_model.nodes;                       % metres
S_mm   = nodes(unique(img.fwd_model.boundary(:)),:) * 1000;
ctr_mm = mean(nodes,1) * 1000;
a_mm   = anchor_xyz_m(:)' * 1000;

dirv = ctr_mm - a_mm;
dirv = dirv / norm(dirv);                           % unit vector, inward

t = (0:0.25:140)';                                  % 0-140 mm along the line
P = a_mm + t.*dirv;

% Distance from each sample point to the nearest scalp node.
d = zeros(size(t));
for i = 1:numel(t)
    d(i) = sqrt(min(sum((S_mm - P(i,:)).^2, 2)));
end

dmax = max(d);

prof = struct('t_mm', t, 'depth_mm', d, 'points_mm', P, ...
              'anchor_mm', a_mm, 'direction', dirv, 'dmax_mm', dmax);
end
