function [pos, names] = a1_make_1010_positions(mesh_name, do_plot)

%A1_MAKE_1010_POSITIONS  Standard EEG 10-10 electrode positions on a Tizzard mesh.

% Aim served: Aim 1 (electrode model) and Aim 2 (electrode configuration).
%
% Follows the geometric construction in EIDORS mk_head_model_adult.m
% (get_elec_pos): slice the scalp with a plane, fit a Fourier curve to the
% outline, and place electrodes at 10% steps along it. The EIDORS nasion and
% nasal-cavity boxes are used, since the Tizzard meshes share its coordinates.
%
%   [pos, names] = a1_make_1010_positions(mesh_name, do_plot)
%     mesh_name  Tizzard mesh to use, default 'SAH262'
%     do_plot    true to plot the positions, default false
%     pos        73 x 3 electrode positions, mm, in the Tizzard frame
%     names      1 x 73 cell array of 10-10 labels
% Imogen Ryan - Honours thesis.

if nargin < 1 || isempty(mesh_name), mesh_name = 'SAH262'; end
if nargin < 2, do_plot = false; end

base = thesis_dir('models','tizzard_meshes');
S = load([base mesh_name '.mat'],'tri','vtx');
fmdl = eidors_obj('fwd_model','tz','nodes',S.vtx,'elems',S.tri);
fmdl.boundary = find_boundary(fmdl);
img = mk_image(fmdl,1);

pos = []; names = {};

%% 1. Mid-sagittal chain, nasion to inion, 10% steps
slc = mdl_slice_mesher(img,[0 inf inf]);
slc.fwd_model.boundary = find_boundary(slc.fwd_model);
b  = unique(slc.fwd_model.boundary);
nd = slc.fwd_model.nodes;
% Exclude the nasal cavity (EIDORS box for this head)
boxA = [-87.6 -63.5; -88.2 -28; 9.6 -26.8; 11.6 -67.5];
b = b(~inpolygon(nd(b,2), nd(b,3), boxA(:,1), boxA(:,2)));
outln = order_loop(nd(b,:), -1);
nasion = [-81.754 -8.157];                       % nasion (mm), EIDORS, this head
pp  = fourier_fit(outln(:,2:3), 51);
% 0.456 of the fitted loop reaches the inion (EIDORS value for this head)
fr  = linspace(0, 0.456, 1001); fr(end) = [];
cap = fourier_fit(pp, fr, nasion);
mid = {'Nz','Fpz','AFz','Fz','FCz','Cz','CPz','Pz','POz','Oz','Iz'};
idxm = [1 100 200 300 400 500 600 700 800 900 1000];
for k = 1:numel(mid)
    [pos,names] = add_e(pos,names,[0 cap(idxm(k),:)],mid{k});
end

%% 2. Circumferential ring through Fpz and Oz
e1 = get_e(pos,names,'Fpz');  e2 = get_e(pos,names,'Oz');
d  = e2 - e1;
y  = e1(2) + (-e1(3)/d(3))*d(2);
z  = e1(3) + (-e1(2)/d(2))*d(3);
slc = mdl_slice_mesher(img,[inf y z]);
slc.fwd_model.boundary = find_boundary(slc.fwd_model);
b = unique(slc.fwd_model.boundary);
outln = order_loop(slc.fwd_model.nodes(b,:), -1);
pp = fourier_fit(outln(:,1:2), 51);
fr = linspace(0,1,21); fr(end) = [];
pl = fourier_fit(pp, fr, e1(:,1:2));
pl(:,3) = e1(3) + d(3)*(pl(:,2) - e1(2))/d(2);
pl([1 11],:) = [];                                % drop Fpz and Oz themselves
ring = {'Fp1','AF7','F7','FT7','T7','TP7','P7','PO7','O1', ...
        'O2','PO8','P8','TP8','T8','FT8','F8','AF8','Fp2'};
for k = 1:numel(ring), [pos,names] = add_e(pos,names,pl(k,:),ring{k}); end

%% 3. Coronal chains, each through a midline and a ring electrode
chains = { ...
  {'AFz','AF7', {'AF3','AF4'},                                     [3 6]}, ...
  {'Fz' ,'F7' , {'F9','F5','F3','F1','F2','F4','F6','F10'},          []}, ...
  {'FCz','FT7', {'FT9','FC5','FC3','FC1','FC2','FC4','FC6','FT10'}, []}, ...
  {'Cz' ,'T7' , {'T9','C5','C3','C1','C2','C4','C6','T10'},          []}, ...
  {'CPz','TP7', {'TP9','CP5','CP3','CP1','CP2','CP4','CP6','TP10'}, []}, ...
  {'Pz' ,'P7' , {'P9','P5','P3','P1','P2','P4','P6','P10'},          []}, ...
  {'POz','PO7', {'PO3','PO4'},                                     [3 6]} };
for c = 1:numel(chains)
    el = coronal_chain(img, pos, names, chains{c}{1}, chains{c}{2});
    if ~isempty(chains{c}{4}), el = el(chains{c}{4},:); end
    nm = chains{c}{3};
    for k = 1:numel(nm), [pos,names] = add_e(pos,names,el(k,:),nm{k}); end
end

%% 4. Report and optional plot
P = fmdl.nodes(unique(fmdl.boundary(:)),:);
dsurf = zeros(size(pos,1),1);
for k = 1:size(pos,1), dsurf(k) = min(sqrt(sum((P-pos(k,:)).^2,2))); end
fprintf('a1_make_1010_positions: %d electrodes on %s\n', size(pos,1), mesh_name);
fprintf('  distance to nearest scalp node: median %.2f mm, max %.2f mm\n', ...
    median(dsurf), max(dsurf));
if max(dsurf) > 6
    warning('a1_make_1010:offSurface', ...
        '%d position(s) more than 6 mm off the scalp -- inspect before use', sum(dsurf>6));
end

if do_plot
    figure('Color','w','Position',[60 60 900 700]);
    trisurf(fmdl.boundary, fmdl.nodes(:,1),fmdl.nodes(:,2),fmdl.nodes(:,3), ...
        'FaceColor',[0.88 0.80 0.74],'EdgeColor','none','FaceAlpha',0.9); hold on
    plot3(pos(:,1),pos(:,2),pos(:,3),'r.','MarkerSize',18);
    for k=1:numel(names)
        text(pos(k,1)*1.03,pos(k,2)*1.03,pos(k,3)*1.03,names{k},'FontSize',7);
    end
    axis equal off; view([-0.4 -1 0.35]); camlight headlight; material dull
    title(sprintf('EEG 10-10 positions on %s (%d electrodes)',mesh_name,size(pos,1)));
end
end

% Helpers
function [pos,names] = add_e(pos,names,p,nm)
pos(end+1,:) = p(:)';  names{end+1} = nm;
end

function e = get_e(pos,names,nm)
i = find(strcmp(names,nm),1);
if isempty(i), error('a1_make_1010:noElec','electrode "%s" not defined yet',nm); end
e = pos(i,:);
end

function el = coronal_chain(img, pos, names, nMid, nLat)
%CORONAL_CHAIN  Electrodes at 10% steps on the arc through two named electrodes.
%  Follows define_plane() in EIDORS mk_head_model_adult.m.
e1 = get_e(pos,names,nMid);  e2 = get_e(pos,names,nLat);
d  = e2 - e1;
y  = e1(2) + (-e1(3)/d(3))*d(2);
z  = e1(3) + (-e1(2)/d(2))*d(3);
slc = mdl_slice_mesher(img,[inf y z]);
slc.fwd_model.boundary = find_boundary(slc.fwd_model);
b  = unique(slc.fwd_model.boundary);
nd = slc.fwd_model.nodes;
% Exclude the nasal cavity (EIDORS x-z box for this head)
boxB = [-16.0484 -19.0245; 22.6613 -15.7946; 23.7903 -92.9910; -20.0806 -91.6990];
b = b(~inpolygon(nd(b,1), nd(b,3), boxB(:,1), boxB(:,2)));
outln = order_loop(nd(b,:), -1);
pp = fourier_fit(outln(:,[1 3]), 51);
fr = linspace(0,1,10001); fr(end) = [];
pl = fourier_fit(pp, fr, e1(:,[1 3]));
% Locate the lateral electrode along the fitted arc
[~,p] = min(sqrt(sum((pl - e2([1 3])).^2, 2)));
p = p/10000;  if p > 0.5, p = 1-p; end
fr = [-5 -3 -2 -1 1 2 3 5]./4 .* p;
el = fourier_fit(pp, fr, e1(:,[1 3]));
% Restore y from z along the e1 to e2 direction, as EIDORS define_plane does.
% Using x instead puts positions well off the scalp.
e = zeros(size(el,1),3);
e(:,[1 3]) = el;
e(:,2) = e1(2) + d(2).*(e(:,3) - e1(3))./d(3);
el = e;
end

function d = thesis_dir(varargin)
% Folder inside this repository
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
