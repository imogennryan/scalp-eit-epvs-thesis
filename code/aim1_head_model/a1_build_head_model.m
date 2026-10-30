function [img, info] = a1_build_head_model(opt)

%A1_BUILD_HEAD_MODEL  Build the four-layer Tizzard head model in SI units.

% Aim served: Aim 1 (model build).
%
% Mesh: Tizzard and Bayford 2007, Physiol. Meas. 28:S163-S182.
% Nodes are converted from mm to m, tissues are found by name, and electrodes
% are built from whole boundary faces.
%
%   [img, info] = a1_build_head_model(opt)
%     opt.mesh                 'SAH262' (default), 'SAH028', 'SAH031', 'SAH035',
%                              'SAH049', 'SAH077' or 'SAH132'
%     opt.sigma_set            'mcdermott_baumann' (default), 'mesh',
%                              'progress_report', 'mcdermott2020' or 'gabriel'
%     opt.z_contact            contact impedance, Ohm.m^2, default 0.079
%                              (1 kOhm on a 10 mm electrode, McDermott 2020)
%     opt.montage              'm19_1020' (default), 'm32_1010' or 'm16_twoplane'
%     opt.n_elec               electrodes on the fallback ring, default 16
%     opt.elec_radius          electrode grouping radius, m, default 8e-3
%                              (used only when elec_target_area_m2 is empty)
%     opt.elec_targets         N x 3 electrode centres, m (overrides the montage)
%     opt.elec_target_area_m2  target electrode area, m^2, default 78.5e-6 (10 mm)
%     opt.amp                  drive current, A, default 100e-6 (McDermott 2020)
%     img                      EIDORS image with one conductivity per element
%     info                     struct of element volumes, tissue indices,
%                              electrode areas and positions, and options used
% Imogen Ryan - Honours thesis.

if nargin < 1, opt = struct(); end
% Defaults
def = struct('mesh','SAH262','sigma_set','mcdermott_baumann','z_contact',0.079, ...
             'n_elec',16,'elec_radius',8e-3,'amp',100e-6, ...
             'elec_target_area_m2',78.5e-6,'elec_targets',[], ...
             'montage','m19_1020');
fn = fieldnames(def);
for k = 1:numel(fn), if ~isfield(opt,fn{k}), opt.(fn{k}) = def.(fn{k}); end, end
% Mesh folder
base = thesis_dir('models','tizzard_meshes');
assert(exist([base 'SAH028.mat'],'file')==2, ...
    'Tizzard meshes not found in %s -- see models/tizzard_meshes/README.md', base);
S = load([base opt.mesh '.mat'],'tri','vtx','materials');
names = strtrim(cellstr(S.materials.names));

% (1) Units: mm to m
nodes = S.vtx / 1000;

fmdl = eidors_obj('fwd_model',['tizzard_' opt.mesh], 'nodes', nodes, 'elems', S.tri);
% Boundary faces (faces on the outside of the head), where electrodes attach
fmdl.boundary = find_boundary(fmdl);
B = fmdl.boundary; Nfb = size(B,1);       % Nfb = number of boundary faces

% (2) Tissue lookup by name
iScalp = region_index(names,'scalp');
iSkull = region_index(names,'skull');
iCSF   = region_index(names,'CSF');
iBrain = region_index(names,'brain');

switch lower(opt.sigma_set)
  case 'mesh'             % values shipped with the mesh, no documented source
    sg = struct('scalp',S.materials.sigma(iScalp),'skull',S.materials.sigma(iSkull), ...
                'csf',  S.materials.sigma(iCSF),  'brain',S.materials.sigma(iBrain));
    src = 'mesh-shipped materials.sigma (no documented source)';
  case 'progress_report'  % IT'IS values used in the progress report
    sg = struct('scalp',0.439,'skull',0.018,'csf',1.88,'brain',0.375);
    src = 'IT''IS via progress report Table 2.2';
  case 'mcdermott2020'    % McDermott 2020, 5 Hz-100 Hz
    sg = struct('scalp',0.23,'skull',0.05,'csf',2.0,'brain',0.10);
    src = 'McDermott 2020 Physiol Meas 41:075010, 5 Hz-100 Hz band';
  case 'mcdermott_baumann'  % thesis set, 100 Hz drive frequency
    % scalp 0.23, skull 0.05, brain 0.10 S/m, McDermott 2020 (5 Hz-100 Hz)
    % CSF 1.79 S/m, Baumann 1997 (human CSF, 37 C, 10 Hz-10 kHz)
    sg = struct('scalp',0.23,'skull',0.05,'csf',1.79,'brain',0.10);
    src = ['McDermott 2020 (5 Hz-100 Hz) with CSF from Baumann 1997 ' ...
           '(human, n=7, 37 C), drive frequency 100 Hz'];
  case 'gabriel'          % proposal set, Gabriel 1996
    sg = struct('scalp',0.44,'skull',0.018,'csf',1.79,'brain',0.10);
    src = 'Proposal, citing Gabriel Lau & Gabriel 1996 (bovine grey matter, in vitro)';
  otherwise, error('unknown sigma_set "%s"', opt.sigma_set);
end

% Electrodes
% Fallback 16-electrode ring, ellipse from EIDORS mk_head_model_adult (m)
RING_Z = 5e-3; ELL_CTR = [0 12e-3]; ELL_A = 80e-3; ELL_B = 95e-3;
% A named montage (default) replaces the fallback ring
if isempty(opt.elec_targets) && ~isempty(opt.montage)
    mf = [thesis_dir('models','montages') 'a1_montages.mat'];
    if exist(mf,'file')
        Mt = load(mf,'M');
        j = find(strcmp({Mt.M.name}, opt.montage), 1);
        if ~isempty(j), opt.elec_targets = Mt.M(j).pos_mm/1000; end
    end
end
if isempty(opt.elec_targets)
    th = linspace(0,2*pi,opt.n_elec+1); th(end) = [];
    targets = [ELL_CTR(1)+ELL_A*cos(th(:)), ELL_CTR(2)+ELL_B*sin(th(:)), ...
               RING_Z*ones(opt.n_elec,1)];
else
    targets = opt.elec_targets;
    opt.n_elec = size(targets,1);
end
% Centroid of every boundary face
fc = squeeze(mean(reshape(fmdl.nodes(B,:), Nfb, 3, 3), 2));

% Stop two electrodes from sharing a face
claimed = false(Nfb,1);
elec = struct('nodes',{},'z_contact',{});
elec_area = zeros(opt.n_elec,1);
[~, face_area_all] = find_electrode_bdy(B, fmdl.nodes, []);   % area of every face
% Build each electrode from whole boundary faces, so EIDORS treats it as a
% complete electrode model with no orphan nodes
for i = 1:opt.n_elec
    % (a) Seed: the unclaimed boundary face nearest this electrode's target
    dd = sqrt(sum((fc - targets(i,:)).^2,2)); dd(claimed) = Inf;
    [~,seedf] = min(dd);

    % (b) Distance of every face from the seed, so the electrode grows as a
    %     compact patch
    d2 = sqrt(sum((fc - fc(seedf,:)).^2,2));
    avail = find(~claimed);
    if isempty(opt.elec_target_area_m2)
        % Fixed grouping radius
        sel = union(find(d2 <= opt.elec_radius & ~claimed), seedf);
    else
        % (c) Fixed target area: add faces in order of distance from the seed
        %     until their total area reaches the target. This keeps electrodes
        %     the same size on coarse and fine meshes.
        [~,ord] = sort(d2(avail));
        cand = avail(ord);                              % faces, nearest first
        cum  = cumsum(face_area_all(cand));             % running total area
        nkeep = find(cum >= opt.elec_target_area_m2, 1);
        if isempty(nkeep), nkeep = numel(cand); end     % ran out of faces
        sel = union(cand(1:nkeep), seedf);
    end
    claimed(sel) = true;
    elec(i).nodes = unique(B(sel,:))';
    elec(i).z_contact = opt.z_contact;
    [~,ba] = find_electrode_bdy(B, fmdl.nodes, elec(i).nodes);
    elec_area(i) = sum(ba);                       % m^2
end
fmdl.electrode = elec;

[stim,msel] = mk_stim_patterns(opt.n_elec,1,'{ad}','{ad}',{'no_meas_current'},opt.amp);
fmdl.stimulation = stim; fmdl.meas_select = msel;
fmdl.solve = @fwd_solve_1st_order; fmdl.system_mat = @system_mat_1st_order;
fmdl.jacobian = @jacobian_adjoint; fmdl.normalize_measurements = 0;
% Ground node: the node nearest the mesh centre, deep in the brain and far
% from every electrode
ctr = mean(fmdl.nodes,1);
[~,fmdl.gnd_node] = min(sum((fmdl.nodes-ctr).^2,2));

% One conductivity per element, set for each tissue by name
img = mk_image(fmdl,1);
img.elem_data(S.materials.ref==iScalp) = sg.scalp;
img.elem_data(S.materials.ref==iSkull) = sg.skull;
img.elem_data(S.materials.ref==iCSF)   = sg.csf;
img.elem_data(S.materials.ref==iBrain) = sg.brain;
img.name = sprintf('%s | sigma: %s', opt.mesh, src);

% Model information returned with the image
a=fmdl.nodes(fmdl.elems(:,1),:); b=fmdl.nodes(fmdl.elems(:,2),:);
c=fmdl.nodes(fmdl.elems(:,3),:); d=fmdl.nodes(fmdl.elems(:,4),:);
info.elem_volume_m3 = abs(dot(b-a,cross(c-a,d-a,2),2))/6;
info.elem_centroid  = (a+b+c+d)/4;
info.ref = S.materials.ref; info.names = names;
info.idx = struct('scalp',iScalp,'skull',iSkull,'csf',iCSF,'brain',iBrain);
info.brain_idx = find(S.materials.ref==iBrain);
info.sigma = sg; info.sigma_source = src;
info.elec_area_m2 = elec_area;
info.elec_centre_m = zeros(opt.n_elec,3);
for i=1:opt.n_elec
    info.elec_centre_m(i,:) = mean(fmdl.nodes(fmdl.electrode(i).nodes,:),1);
end
info.elec_target_m = targets;
info.elec_snap_err_mm = sqrt(sum((info.elec_centre_m-targets).^2,2))*1000;
info.opt = opt;
info.units = 'nodes in m; sigma in S/m; z_contact in Ohm.m^2; current in A; voltages in V';
end

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
