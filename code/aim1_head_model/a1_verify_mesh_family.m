%A1_VERIFY_MESH_FAMILY  Check all seven Tizzard adult meshes before using them.

% Aim served: Aim 1 (model build and validation).
%
% Checks element and node counts against the published EIDORS table, the tissue
% order in each file, inverted elements, tissue volumes and element quality.
% SAH031 stores its tissues in a different order from the other six meshes, so
% tissues are always found by name with region_index, never by number.
%
% Inputs : SAH028/031/035/049/077/132/262 .mat
% Outputs: results/aim1/a1_mesh_family_verification.csv and .mat
% Imogen Ryan - Honours thesis.

clear; close all force;
addpath(fileparts(mfilename('fullpath')));
assert(exist('calc_mesh_quality','file')==2, 'Start EIDORS first (needed for calc_mesh_quality)');
% Mesh folder
base = thesis_dir('models','tizzard_meshes');
assert(exist([base 'SAH028.mat'],'file')==2, ...
    'Tizzard meshes not found in %s -- see models/tizzard_meshes/README.md', base);
outd = thesis_dir('results','aim1');
if ~exist(outd,'dir'), mkdir(outd); end

names_pub = {'SAH028','SAH031','SAH035','SAH049','SAH077','SAH132','SAH262'};
elems_pub = [27684 31111 35073 49172 77196 131672 262368];   % EIDORS data_contrib page
nodes_pub = [ 5226  5939  6556  9142 14188  24688  48394];   % EIDORS data_contrib page

fprintf('\n=========== A1: TIZZARD MESH FAMILY VERIFICATION ===========\n\n');
T = table(); Tout = table(); allok = true; ref_order = {};

for m = 1:numel(names_pub)
    f = [base names_pub{m} '.mat'];
    if ~exist(f,'file')
        fprintf('%-8s MISSING\n', names_pub{m}); allok=false; continue
    end
    S = load(f,'tri','vtx','materials');
    nm = strtrim(cellstr(S.materials.names));
    ok_elem = size(S.tri,1)  == elems_pub(m);
    ok_node = size(S.vtx,1)  == nodes_pub(m);

    % Order string, e.g. 'skull|brain|scalp|CSF'
    ordstr = strjoin(nm(:)','|');
    if m==1, ref_order = ordstr; end
    ok_order = strcmpi(ordstr, ref_order);

    % Per-region element counts must sum to the total
    cnt = arrayfun(@(k) sum(S.materials.ref==k), 1:numel(nm));
    ok_sum = sum(cnt) == size(S.tri,1);

    % Look up tissues by name, since SAH031 orders its regions differently
    iBrain = region_index(nm,'brain');
    iSkull = region_index(nm,'skull');
    iScalp = region_index(nm,'scalp');
    iCSF   = region_index(nm,'CSF');

    % Element quality: a tetrahedron with negative or zero signed volume is
    % inverted or flat, which would make the stiffness matrix wrong.
    a=S.vtx(S.tri(:,1),:); b=S.vtx(S.tri(:,2),:);
    c=S.vtx(S.tri(:,3),:); d=S.vtx(S.tri(:,4),:);
    sv  = dot(b-a, cross(c-a,d-a,2), 2)/6;      % signed volume (mm^3)
    vol = abs(sv);
    n_neg = sum(sv<=0);                          % inverted or flat elements

    % Conductivities shipped with each mesh
    sig = S.materials.sigma(:)';

    ok = ok_elem && ok_node && ok_order && ok_sum;
    allok = allok && ok;
    fprintf('%-8s elems %7d (%s)  nodes %6d (%s)  order %-28s %s\n', ...
        names_pub{m}, size(S.tri,1), tf(ok_elem), size(S.vtx,1), tf(ok_node), ...
        ordstr, tf(ok_order));
    fprintf('         regions sum %s | inverted tets %d | sigma [%s]\n', ...
        tf(ok_sum), n_neg, num2str(sig,'%.3f '));
    vBrain = vol(S.materials.ref==iBrain);
    vSkull = vol(S.materials.ref==iSkull);
    fprintf('         idx brain=%d skull=%d scalp=%d CSF=%d | sigma_brain %.3f sigma_skull %.3f\n', ...
        iBrain, iSkull, iScalp, iCSF, sig(iBrain), sig(iSkull));
    fprintf('         brain median elem %8.2f mm^3 | skull median %8.2f mm^3 | brain/skull %.2f\n', ...
        median(vBrain), median(vSkull), median(vBrain)/median(vSkull));

    T = [T; table(string(names_pub{m}), size(S.tri,1), size(S.vtx,1), ...
        ok_elem, ok_node, ok_order, n_neg, iBrain, iSkull, ...
        median(vBrain), median(vSkull), sum(vBrain)/1e6, ...
        'VariableNames',{'mesh','n_elems','n_nodes','count_ok','node_ok', ...
        'order_matches_SAH028','inverted_tets','idx_brain','idx_skull', ...
        'brain_med_mm3','skull_med_mm3','brain_vol_L'})];

    % Table for the CSV: published vs counted, tissue index and volume per layer
    tissue_mm3 = @(i) sum(vol(S.materials.ref==i));     % mm^3

    % Mesh quality (element shape), EIDORS calc_mesh_quality. eta scores each
    % element from 0 (flat) to 1 (regular tetrahedron)
    qmdl = eidors_obj('fwd_model','q','nodes',S.vtx,'elems',S.tri, ...
                      'boundary',find_boundary(S.tri));
    Q = calc_mesh_quality(qmdl);
    eta = Q.tet.eta;
    label = string(names_pub{m});
    if strcmp(names_pub{m},'SAH262'), label = label + " (used)"; end   % working mesh
    Tout = [Tout; table(label, ...
        elems_pub(m), size(S.tri,1), yesno(ok_elem), ...
        nodes_pub(m), size(S.vtx,1), yesno(ok_node), ...
        yesno(ok_sum), n_neg, string(ordstr), ...
        iScalp, iSkull, iCSF, iBrain, ...
        tissue_mm3(iScalp), tissue_mm3(iSkull), tissue_mm3(iCSF), tissue_mm3(iBrain), ...
        sum(vol), min(vol), median(vol), max(vol), ...
        median(vBrain), median(vSkull), median(vBrain)/median(vSkull), ...
        median(eta), min(eta), sum(eta<0.1), sum(eta<0.3), ...
        'VariableNames',{'mesh', ...
        'published_elements','counted_elements','elements_match', ...
        'published_nodes','counted_nodes','nodes_match', ...
        'all_elements_labelled','inverted_elements','tissue_order', ...
        'idx_scalp','idx_skull','idx_CSF','idx_brain', ...
        'scalp_mm3','skull_mm3','CSF_mm3','brain_mm3','head_mm3', ...
        'elem_vol_min_mm3','elem_vol_median_mm3','elem_vol_max_mm3', ...
        'brain_elem_median_mm3','skull_elem_median_mm3','brain_to_skull_elem_ratio', ...
        'quality_median','quality_worst','quality_below_0p1','quality_below_0p3'})];
end

% Tissue volumes compared with the finest mesh (SAH262), in percent
for tissue = {'scalp','skull','CSF','brain','head'}
    v = Tout.([tissue{1} '_mm3']);
    Tout.([tissue{1} '_diff_from_SAH262_pct']) = 100*(v - v(end))/v(end);
end

fprintf('\n----------------------------------------------------------\n');
fprintf('Element and node counts vs the published table: %s\n', ...
    tf(all(T.count_ok) && all(T.node_ok)));
fprintf('Inverted tetrahedra anywhere in the family      : %d\n', sum(T.inverted_tets));
if all(T.order_matches_SAH028)
    fprintf('Region ordering identical across all meshes     : OK\n');
else
    bad = T.mesh(~T.order_matches_SAH028);
    fprintf('Region ordering DIFFERS in: %s\n', strjoin(cellstr(bad)', ', '));
    fprintf('  -> This is NOT corruption. Each file pairs names with sigma correctly.\n');
    fprintf('  -> It DOES mean tissue indices are mesh-specific.\n');
    fprintf('  -> ALWAYS resolve tissues by NAME. Never hard-code ref==2.\n');
end
fprintf('\nBrain-to-skull element volume ratio by mesh (Tizzard 2005 states the\n');
fprintf('brain should be meshed AS FINE AS the skull, i.e. ratio near 1 or below):\n');
for m=1:height(T)
    fprintf('   %-8s %7d elems   ratio %.2f\n', T.mesh(m), T.n_elems(m), ...
        T.brain_med_mm3(m)/T.skull_med_mm3(m));
end
writetable(Tout,[outd 'a1_mesh_family_verification.csv']);
save([outd 'a1_mesh_family_verification.mat'],'T','Tout','ref_order','allok');
fprintf('SAVED %s\n', [outd 'a1_mesh_family_verification.csv']);

function s = tf(x); if x, s='OK'; else, s='**DIFFERS**'; end; end
function s = yesno(x); if x, s="TRUE"; else, s="FALSE"; end; end

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
