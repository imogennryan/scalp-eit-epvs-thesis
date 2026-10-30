function [f, place] = a2_place_perturbation(img, info, opt)

%A2_PLACE_PERTURBATION  Place one simulated accumulation in the brain at a set depth below a scalp point.

% Aim served: Aim 2.
%
% The sphere is placed on a straight line from a scalp point toward the head
% centre, so depth changes while lateral position stays fixed. Depth is the
% distance from the sphere centre to the nearest scalp node. The sphere is
% mapped onto the mesh as the fraction of each element inside it. This
% function returns fractions only; conductivity is set by a2_mix_conductivity_L.
%
%   [f, place] = a2_place_perturbation(img, info, opt)
%     opt.mode          'radial' (default), or 'nearest' to show the depth confound
%     opt.depth_mm      target depth below the scalp, mm
%     opt.diam_mm       sphere diameter, mm (default 10)
%     opt.anchor_xyz_m  1x3 scalp point the line starts from, m
%     opt.anchor_elec   electrode index to start from, used only when
%                       anchor_xyz_m is empty (default 10). Indices differ
%                       between montages, so use anchor_xyz_m across montages.
%     opt.brain_only    true (default), perturb brain tissue only
%     f                 n_elems x 1, fraction of each element inside the sphere, 0 to 1
%     place             struct recording where the sphere was placed
% Imogen Ryan - Honours thesis.

if nargin < 3, opt = struct(); end
def = struct('mode','radial','depth_mm',25,'diam_mm',10, ...
             'anchor_elec',10,'anchor_xyz_m',[],'brain_only',true);
fn = fieldnames(def);
for k = 1:numel(fn), if ~isfield(opt,fn{k}), opt.(fn{k}) = def.(fn{k}); end, end

fmdl  = img.fwd_model;
nodes = fmdl.nodes;                       % node positions (m)
B     = fmdl.boundary;                    % scalp surface faces
scalp_nodes = unique(B(:));               % every node on the outer surface
S_mm  = nodes(scalp_nodes,:) * 1000;      % scalp nodes in mm
cen_mm = info.elem_centroid * 1000;       % element centroids in mm

% Choose the sphere centre
switch lower(opt.mode)

    case 'radial'
        % Straight line from the scalp point toward the head centre. Moving
        % along it changes depth while lateral position stays fixed. A
        % physical scalp point keeps the same line across montages.
        if ~isempty(opt.anchor_xyz_m)
            anchor_mm = opt.anchor_xyz_m(:)' * 1000;
        else
            assert(opt.anchor_elec <= numel(fmdl.electrode), ...
                ['anchor_elec %d does not exist in this model (%d electrodes). ' ...
                 'Electrode indices are montage-specific -- pass anchor_xyz_m instead.'], ...
                opt.anchor_elec, numel(fmdl.electrode));
            anchor_mm = info.elec_centre_m(opt.anchor_elec,:) * 1000;
        end
        ctr_mm    = mean(nodes,1) * 1000;                 % head centroid
        dirv      = ctr_mm - anchor_mm;
        dirv      = dirv / norm(dirv);                    % unit vector, inward

        % Step inward until the point first reaches the target depth below
        % the scalp. The line passes through the head, so each depth is
        % reached twice. The first crossing keeps the region on the near side.
        t  = (0:0.25:120)';        % 0 to 120 mm along the line (a2_reachable_depth scans 0 to 140 mm)
        P  = anchor_mm + t.*dirv;
        dP = depth_of(P, S_mm);                           % depth of each point (mm)
        it = find(dP >= opt.depth_mm, 1, 'first');        % first crossing only
        if isempty(it)
            % Name the anchor that was used in the error message
            if ~isempty(opt.anchor_xyz_m)
                whichanchor = sprintf('scalp point [%.1f %.1f %.1f] mm', ...
                                      opt.anchor_xyz_m(:)'*1000);
            else
                whichanchor = sprintf('electrode index %d', opt.anchor_elec);
            end
            error('a2_place_perturbation:tooDeep', ...
                ['Depth %.1f mm is never reached along the radial line from %s ' ...
                 '(maximum reachable is %.1f mm). Head shape limits how deep a ' ...
                 'line from this position can go; see a2_reachable_depth.'], ...
                opt.depth_mm, whichanchor, max(dP));
        end
        centre_mm = P(it,:);

    case 'nearest'
        % Brain element nearest the target depth, at any lateral position.
        % Kept to show the depth confound. Not used for results.
        dEl = depth_of(cen_mm(info.brain_idx,:), S_mm);
        [~, ib] = min(abs(dEl - opt.depth_mm));
        centre_mm = cen_mm(info.brain_idx(ib),:);

    otherwise
        error('a2_place_perturbation:mode','unknown mode "%s"', opt.mode);
end

% Map the sphere onto the mesh as element fractions
% mk_c2f_circ_mapping takes [x; y; z; radius] in metres. Diameter is halved for the radius.
xyzr = [centre_mm(:)/1000; (opt.diam_mm/2)/1000];
map  = mk_c2f_circ_mapping(fmdl, xyzr);
f    = full(map(:,1));
f(~isfinite(f)) = 0;                       % guard against NaN
f    = max(min(f,1), 0);                   % clamp to a physical fraction

% Keep the perturbation inside brain tissue
if opt.brain_only
    inbrain = false(size(f));
    inbrain(info.brain_idx) = true;
    f(~inbrain) = 0;
end

% Record what was placed
place = struct();
place.mode              = opt.mode;
place.centre_mm         = centre_mm;
place.depth_mm_target   = opt.depth_mm;
place.depth_mm_actual   = depth_of(centre_mm, S_mm);   % recomputed, not assumed
place.diam_mm           = opt.diam_mm;
place.vol_mm3           = sum(f .* info.elem_volume_m3) * 1e9;
place.vol_mm3_nominal   = 4/3*pi*(opt.diam_mm/2)^3;
place.n_elems_touched   = sum(f > 0);
place.max_fraction      = max(f);
% Share of the intended sphere that landed in brain tissue. Near the surface
% part of the sphere falls in scalp, skull or CSF and is removed, so a sweep
% can reject placements that do not fit.
place.vol_fraction_delivered = place.vol_mm3 / place.vol_mm3_nominal;
place.fits_in_brain          = place.vol_fraction_delivered > 0.99;
% Distance to the nearest electrode, recorded for every placement
d_elec = sqrt(sum((info.elec_centre_m*1000 - centre_mm).^2, 2));
place.dist_to_nearest_elec_mm = min(d_elec);
place.anchor_elec       = opt.anchor_elec;
place.anchor_mm         = anchor_mm_out(opt, info);

end

function a = anchor_mm_out(opt, info)
% The scalp point the line started from (mm).
if ~isempty(opt.anchor_xyz_m), a = opt.anchor_xyz_m(:)' * 1000;
else,                          a = info.elec_centre_m(opt.anchor_elec,:) * 1000;
end
end

function d = depth_of(P, S_mm)
% Distance from each point in P to the nearest scalp node (mm).
d = zeros(size(P,1),1);
for i = 1:size(P,1)
    d(i) = sqrt(min(sum((S_mm - P(i,:)).^2, 2)));
end
end
