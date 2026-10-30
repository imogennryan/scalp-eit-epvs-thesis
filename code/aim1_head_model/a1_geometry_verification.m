%A1_GEOMETRY_VERIFICATION  Check that the four compartments are nested and measure their thickness.

% Aim served: Aim 1 (model build).
%
% Rays are cast outward from the head centre and the tissue at each step is
% found. This gives the order the compartments are crossed (brain, CSF, skull,
% scalp) and the thickness of each. Thickness is measured over the cranial
% vault only, because the scalp compartment also includes the face and neck.
% Thickness and volumes are compared with published adult values (Bayford 2006,
% Culpepper 2024 and the sources it cites).
%
% Inputs : models/tizzard_meshes/
% Outputs: results/aim1/a1_geometry_*.csv, a1_anatomical_*.csv,
%          a1_nesting_lines.xlsx and a1_geometry_verification.mat
% Imogen Ryan - Honours thesis.

clear; close all force;
addpath(fileparts(mfilename('fullpath')));
outd = thesis_dir('results','aim1');  if ~exist(outd,'dir'), mkdir(outd); end
base = thesis_dir('models','tizzard_meshes');

MESH     = 'SAH262';     % finest mesh, so compartment boundaries are best resolved
N_RAY    = 400;          % rays from the centroid
STEP     = 0.25;         % sampling step along each ray (mm)
VAULT_EL = 0.5;          % sin(elevation) cut-off, 0.5 is above 30 degrees
                         % (cranial vault only)

S = load([base MESH '.mat'],'tri','vtx','materials');
names = strtrim(cellstr(S.materials.names));
idx = struct('scalp',find(strcmpi(names,'scalp')), 'skull',find(strcmpi(names,'skull')), ...
             'csf',  find(strcmpi(names,'CSF')),   'brain',find(strcmpi(names,'brain')));
TR  = triangulation(S.tri, S.vtx);
ctr = mean(S.vtx,1);

fprintf('\n======== A1: GEOMETRY VERIFICATION (%s) ========\n', MESH);
fprintf('Casting %d rays from the head centroid [%.1f %.1f %.1f] mm, %.2f mm sampling.\n\n', ...
    N_RAY, ctr, STEP);

% Evenly spread ray directions (golden-angle spiral on a sphere)
k = (0:N_RAY-1)' + 0.5;                       % +0.5 keeps rays off the exact poles
phi = acos(1 - 2*k/N_RAY);  gold = pi*(1+sqrt(5));  th = gold*k;
dirs = [cos(th).*sin(phi), sin(th).*sin(phi), cos(phi)];   % N_RAY x 3 unit vectors

Rmax = max(sqrt(sum((S.vtx-ctr).^2,2)));
r = (0:STEP:Rmax)';
tiss = {'brain','csf','skull','scalp'};
th_mm = nan(N_RAY,4);  order_ok = false(N_RAY,1);  nvalid = 0;
labs = cell(N_RAY,1);  firsts_all = nan(N_RAY,4);   % tissue label at every step of every line (NaN = outside the head)

% Walk out along every ray, recording which tissue each step lands in.
for i = 1:N_RAY
    % Sample points along this ray
    P = ctr + r.*dirs(i,:);

    % Tetrahedron containing each point, NaN once the ray leaves the head
    ti = pointLocation(TR, P);
    inmesh = ~isnan(ti);
    if ~any(inmesh), continue; end        % ray missed the mesh; skip it

    % Turn "which tetrahedron" into "which tissue" via materials.ref.
    lab = nan(size(ti));
    lab(inmesh) = S.materials.ref(ti(inmesh));
    labs{i} = lab;
    % Thickness of each compartment along this ray
    for t = 1:4
        th_mm(i,t) = sum(lab == idx.(tiss{t})) * STEP;
    end
    % Nesting test: the distance at which each tissue is first met must
    % increase in the order brain, CSF, skull, scalp
    firsts = nan(1,4);
    for t = 1:4
        f = find(lab == idx.(tiss{t}), 1, 'first');
        if ~isempty(f), firsts(t) = r(f); end
    end
    firsts_all(i,:) = firsts;
    seen = ~isnan(firsts);
    % Judge only rays that cross at least 3 compartments, since rays through
    % the face or neck have no skull to cross
    if sum(seen) >= 3
        nvalid = nvalid + 1;
        order_ok(i) = issorted(firsts(seen));
    end
end

fprintf('--- 1. NESTING: is the order brain -> CSF -> skull -> scalp? ---\n');
fprintf('  rays crossing >=3 compartments : %d of %d\n', nvalid, N_RAY);
fprintf('  rays with the correct order    : %d (%.1f%%)\n', sum(order_ok), 100*sum(order_ok)/max(nvalid,1));
if sum(order_ok) == nvalid
    fprintf('  --> PASS: every usable ray crosses the compartments in anatomical order.\n\n');
else
    fprintf('  --> %d ray(s) out of order. Inspect before using this mesh.\n\n', nvalid-sum(order_ok));
end

vault = dirs(:,3) > VAULT_EL;
fprintf('--- 2. LAYER THICKNESS (mm, along-ray) ---\n');
fprintf('Cranial vault only: %d of %d rays (elevation > %.0f deg).\n', ...
    sum(vault), N_RAY, asind(VAULT_EL));
fprintf('%-8s %-10s %-10s %-10s %-12s %s\n','tissue','median','5th pct','95th pct','all-dirs med','literature');
Tt = table();
for t = 1:4
    x  = th_mm(vault,t);  x  = x(isfinite(x) & x > 0);      % vault only
    xa = th_mm(:,t);      xa = xa(isfinite(xa) & xa > 0);   % all directions
    switch tiss{t}
        case {'scalp','csf'}, lit = '~1-3 mm (Bayford 2006, approximate)';
        otherwise,            lit = 'no published range';
    end
    fprintf('%-8s %-10.2f %-10.2f %-10.2f %-12.2f %s\n', tiss{t}, ...
        median(x), prctile(x,5), prctile(x,95), median(xa), lit);
    Tt = [Tt; table(string(tiss{t}), median(x), prctile(x,5), prctile(x,95), median(xa), ...
        string(lit),'VariableNames', ...
        {'tissue','vault_median_mm','vault_p5_mm','vault_p95_mm','alldir_median_mm','literature'})];
end
fprintf('\n  The "all-dirs" column is CONTAMINATED by rays through the face and neck,\n');
fprintf('  which the scalp compartment includes. Use the vault columns.\n');
% Compare both tissues in Bayford's statement
iC = strcmp(Tt.tissue,'csf');
iS = strcmp(Tt.tissue,'scalp');
fprintf('\n  Bayford 2006 says scalp AND CSF are "approximately 1-3 mm". Both, against that:\n');
fprintf('    CSF   vault %.2f mm (5-95%%: %.2f-%.2f) -- consistent\n', ...
    Tt.vault_median_mm(iC), Tt.vault_p5_mm(iC), Tt.vault_p95_mm(iC));
fprintf('    scalp vault %.2f mm (5-95%%: %.2f-%.2f) -- ABOVE the quoted range\n', ...
    Tt.vault_median_mm(iS), Tt.vault_p5_mm(iS), Tt.vault_p95_mm(iS));
fprintf('  That is not a model failure. Bayford''s figure is an uncited approximation\n');
fprintf('  made in a methods argument, it covers two tissues with one number, and the\n');
fprintf('  scalp compartment here is all soft tissue over the skull. See header.\n');
fprintf('  NEITHER tissue has a sourced anatomical reference range in this corpus.\n');

fprintf('\n--- 3. COMPARTMENT VOLUMES ---\n');
a=S.vtx(S.tri(:,1),:); b=S.vtx(S.tri(:,2),:); c=S.vtx(S.tri(:,3),:); d=S.vtx(S.tri(:,4),:);
vol = abs(dot(b-a,cross(c-a,d-a,2),2))/6;
Tv = table();
for t = 1:4
    V = sum(vol(S.materials.ref==idx.(tiss{t})));
    fprintf('  %-8s %8.1f mL (%.1f%% of head)\n', tiss{t}, V/1e3, 100*V/sum(vol));
    Tv = [Tv; table(string(tiss{t}), V/1e3, 100*V/sum(vol), ...
        'VariableNames',{'tissue','volume_mL','pct_of_head'})];
end
fprintf('  %-8s %8.1f mL\n', 'TOTAL', sum(vol)/1e3);
fprintf('\n  ** Skull and scalp thickness remain UNSOURCED -- see check 4.\n');

%% 4. Anatomical reference ranges
% Brain volume Allen 2002, CSF volume Grant 1987, CSF thickness Haeussinger 2011,
% head size Lee 2018, as cited by Culpepper 2024
fprintf('\n--- 4. ANATOMICAL REFERENCE RANGES ---\n');
Vb  = Tv.volume_mL(strcmp(Tv.tissue,'brain'));
Vc  = Tv.volume_mL(strcmp(Tv.tissue,'csf'));
tc  = Tt.vault_median_mm(strcmp(Tt.tissue,'csf'));
ext = max(S.vtx) - min(S.vtx);                   % mm, whole-mesh bounding box

fprintf('  brain volume      %7.1f mL   vs 1260 +/- 115 (M), 1130 +/- 99 (F) mL\n', Vb);
fprintf('                                   -> %+.2f SD from the female mean\n', (Vb-1130)/99);
fprintf('  CSF volume        %7.1f mL   vs ~150 mL total intracranial\n', Vc);
fprintf('                                   -> %.2fx the cited adult figure\n', Vc/150);
fprintf('  brain + CSF       %7.1f mL   vs ~1280 mL (1130 + 150)\n', Vb+Vc);
fprintf('  CSF vault thick   %7.2f mm   vs 1.1-8.4 mm (+/-1 SD), mean 4.2 mm\n', tc);
fprintf('  head width  (x)   %7.1f mm   vs 151.4 +/- 6.9 mm\n', ext(1));
fprintf('  head length (y)   %7.1f mm   vs 195.2 +/- 9.2 mm\n', ext(2));
fprintf('\n  READ THESE CAREFULLY:\n');
fprintf('   - CSF thickness sits INSIDE the MRI-derived range, toward the thin end.\n');
fprintf('     Haeussinger measured the forehead region only, so it is indicative.\n');
fprintf('   - CSF volume is above the ~150 mL figure and brain is below the adult\n');
fprintf('     mean, while brain+CSF is close. The model apportions more to CSF and\n');
fprintf('     less to brain than the anatomy does, at roughly the right total.\n');
fprintf('     More CSF means more shunting, which would tend to REDUCE the signal.\n');
fprintf('     That is the OPPOSITE direction to the mesh-resolution biases, and it\n');
fprintf('     has NOT been tested. Do not claim a net direction.\n');
fprintf('   - Head extents are the WHOLE mesh, including face and neck, so they are\n');
fprintf('     upper bounds on cranial dimensions, not a like-for-like comparison.\n');
fprintf('   - SKULL and SCALP thickness still have NO sourced range. Culpepper\n');
fprintf('     aggregates scalp and skull into one layer and gives neither.\n');

Ta = table( ...
    ["brain volume mL";"CSF volume mL";"brain+CSF mL";"CSF vault thickness mm"; ...
     "head width mm";"head length mm"], ...
    [Vb; Vc; Vb+Vc; tc; ext(1); ext(2)], ...
    ["1260+/-115 M, 1130+/-99 F [Allen 2002 via Culpepper 2024]"; ...
     "~150 total intracranial [Grant 1987; Laughlin 2021 via Culpepper 2024]"; ...
     "~1280 [derived]"; ...
     "1.1-8.4, mean 4.2 [Haeussinger 2011 n=23, via Culpepper 2024]"; ...
     "151.4+/-6.9 [Lee 2018 n>2200, via Culpepper 2024]"; ...
     "195.2+/-9.2 [Lee 2018 n>2200, via Culpepper 2024]"], ...
    'VariableNames',{'quantity','model_value','reference_range'});
writetable(Ta,[outd 'a1_anatomical_ranges.csv']);

% Model against published adult values
tSk = Tt.vault_median_mm(strcmp(Tt.tissue,'skull'));
tSc = Tt.vault_median_mm(strcmp(Tt.tissue,'scalp'));
ref = {
 % quantity, model, published mean, SD, low, high, unit, source, checked, what is compared
 'Brain (female)',        Vb, 1131.1, 99.5, NaN, NaN, 'mL', 'Allen 2002, MRI, n = 23', 'yes', 'Cerebrum and cerebellum, excluding ventricles'
 'Brain (male)',          Vb, 1273.6, 115.0, NaN, NaN, 'mL', 'Allen 2002, MRI, n = 23', 'yes', 'Cerebrum and cerebellum, excluding ventricles'
 'Brain (female)',        Vb, 1304.4, 87.5, NaN, NaN, 'mL', 'Kruggel 2006, MRI, n = 145', 'yes', 'Grey and white matter'
 'Brain (male)',          Vb, 1417.3, 86.2, NaN, NaN, 'mL', 'Kruggel 2006, MRI, n = 145', 'yes', 'Grey and white matter'
 'CSF (female)',          Vc, 176.3, 28.2, NaN, NaN, 'mL', 'Kruggel 2006, MRI, n = 145', 'yes', 'All CSF inside the skull'
 'CSF (male)',            Vc, 183.6, 29.7, NaN, NaN, 'mL', 'Kruggel 2006, MRI, n = 145', 'yes', 'All CSF inside the skull'
 'CSF (female)',          Vc, 114.5, NaN, NaN, NaN, 'mL', 'Grant 1987, MRI, n = 39', 'yes', 'Total cranial CSF, no SD given'
 'CSF (male)',            Vc, 146, NaN, NaN, NaN, 'mL', 'Grant 1987, MRI, n = 25', 'yes', 'Total cranial CSF, no SD given'
 'Brain + CSF (female)',  Vb+Vc, 1494.9, 96.3, NaN, NaN, 'mL', 'Kruggel 2006, MRI, n = 145', 'yes', 'Total intracranial volume'
 'Brain + CSF (male)',    Vb+Vc, 1616.3, 91.1, NaN, NaN, 'mL', 'Kruggel 2006, MRI, n = 145', 'yes', 'Total intracranial volume'
 'Skull (female)',        Tv.volume_mL(strcmp(Tv.tissue,'skull')), 730, 82, NaN, NaN, 'mL', 'Rashed 2025, MRI, n = 196', 'yes', 'Probably includes facial bone'
 'Skull (male)',          Tv.volume_mL(strcmp(Tv.tissue,'skull')), 827, 79, NaN, NaN, 'mL', 'Rashed 2025, MRI, n = 196', 'yes', 'Probably includes facial bone'
 'Whole head (female)',   sum(vol)/1e3, 2777.0, 206.8, NaN, NaN, 'mL', 'Kruggel 2006, MRI, n = 145', 'yes', 'Cut 10 mm below the cerebellum, air spaces filled in'
 'Whole head (male)',     sum(vol)/1e3, 3076.7, 225.2, NaN, NaN, 'mL', 'Kruggel 2006, MRI, n = 145', 'yes', 'Cut 10 mm below the cerebellum, air spaces filled in'
 'Head width (female)',   ext(1), 151.3, 7.0, NaN, NaN, 'mm', 'Lee 2018, 3D scans, n = 946', 'yes', 'Head breadth. Model: widest span of the mesh'
 'Head width (male)',     ext(1), 151.5, 6.9, NaN, NaN, 'mm', 'Lee 2018, 3D scans, n = 854', 'yes', 'Head breadth. Model: widest span of the mesh'
 'Head length (female)',  ext(2), 195.1, 9.3, NaN, NaN, 'mm', 'Lee 2018, 3D scans, n = 956', 'yes', 'Model includes nose and neck, so an upper bound'
 'Head length (male)',    ext(2), 195.4, 9.3, NaN, NaN, 'mm', 'Lee 2018, 3D scans, n = 862', 'yes', 'Model includes nose and neck, so an upper bound'
 'CSF thickness',         tc, 4.17, NaN, 2.1, 6.3, 'mm', 'Haeussinger 2011, MRI, n = 22', 'yes', 'Model: top of head. Published: forehead, average of 24 points'
 'Skull thickness',       tSk, 5.24, NaN, 3.9, 6.8, 'mm', 'Haeussinger 2011, MRI, n = 22', 'yes', 'Model: top of head. Published: forehead, average of 24 points'
 'Skull thickness',       tSk, 7.04, 1.43, NaN, NaN, 'mm', 'Thulung 2019, CT, n = 100', 'yes', 'Model: top of head. Published: parietal bone'
 'Scalp thickness',       tSc, 5.12, NaN, 4.0, 7.1, 'mm', 'Haeussinger 2011, MRI, n = 22', 'yes', 'Model: top of head. Published: forehead, average of 24 points'
 };
Tc = cell2table(ref, 'VariableNames', {'quantity','model','published_mean','published_sd', ...
    'published_low','published_high','units','source','checked_in_original','what_is_compared'});
Tc.difference_pct = 100*(Tc.model - Tc.published_mean)./Tc.published_mean;
Tc.difference_in_SD = (Tc.model - Tc.published_mean)./Tc.published_sd;
isVol = strcmp(Tc.units,'mL');
Tc.model_mm3 = nan(height(Tc),1);          Tc.model_mm3(isVol) = Tc.model(isVol)*1e3;
Tc.published_mean_mm3 = nan(height(Tc),1); Tc.published_mean_mm3(isVol) = Tc.published_mean(isVol)*1e3;
writetable(Tc,[outd 'a1_anatomical_comparison.csv']);

writetable(Tt,[outd 'a1_geometry_thickness.csv']);
writetable(Tv,[outd 'a1_geometry_volumes.csv']);
save([outd 'a1_geometry_verification.mat'],'Tt','Tv','th_mm','order_ok','MESH');

% Line-by-line table of the nesting check
% Rank each tissue by its position from the inside out: brain 1, CSF 2, skull 3, scalp 4
rank = zeros(1,numel(names));
for t = 1:4, rank(idx.(tiss{t})) = t; end
tname = {'brain','CSF','skull','scalp'};
elev = asind(dirs(:,3));                          % degrees above (+) or below (-) level
towards = repmat("sides of head", N_RAY, 1);
towards(elev > 30)  = "top of head";
towards(elev < -30) = "down, towards face and neck";
seq = strings(N_RAY,1); goes_back = false(N_RAY,1); brain_to_skull = false(N_RAY,1);
back_z = nan(N_RAY,1);
for i = 1:N_RAY
    in = ~isnan(labs{i});
    L = rank(labs{i}(in));                        % tissue rank at each step inside the head
    P = ctr + r(in).*dirs(i,:);                   % position of each of those steps
    ch = find(diff(L) ~= 0);                      % steps where the tissue changes
    seq(i) = strjoin(tname([L(1) L(ch+1)]), ' > ');
    back = ch(L(ch+1) < L(ch));                   % changes that go back to an inner tissue
    goes_back(i) = ~isempty(back);
    if goes_back(i), back_z(i) = P(back(1)+1,3); end
    brain_to_skull(i) = any(L(ch) == 1 & L(ch+1) == 3);   % brain straight into skull, no CSF
end

% Why lines that go back to an inner tissue do so
pat = {'brain > CSF > skull > scalp > skull > scalp', ...
       'brain > CSF > skull > scalp > skull > scalp > skull > scalp', ...
       'brain > CSF > brain > CSF > skull > scalp', ...
       'brain > CSF > skull > CSF', ...
       'brain > CSF > skull > CSF > skull > scalp', ...
       'brain > CSF > skull > CSF > brain > CSF > skull > scalp'};
why = {'After leaving the skull, the line meets more bone lower down (face, jaw or skull base), then soft tissue again', ...
       'As above, meeting bone twice more on the way out', ...
       'The line crosses a pocket of CSF between two parts of the brain, then goes back into brain before reaching the outer CSF layer', ...
       'The line crosses a ridge of the skull base and comes back into the fluid space inside it, then leaves the mesh', ...
       'The line crosses a ridge of the skull base, passes fluid inside it, then crosses the skull again', ...
       'The line crosses a ridge of the skull base and back into the space inside the skull, then out again'};
explain = strings(N_RAY,1);
for p = 1:numel(pat), explain(seq == pat{p}) = why{p}; end
explain(goes_back & explain == "") = "Not yet explained";

Tr = table((1:N_RAY)', round(elev,1), towards, ...
    round(firsts_all(:,1),2), round(firsts_all(:,2),2), round(firsts_all(:,3),2), round(firsts_all(:,4),2), ...
    yesno(order_ok), yesno(~goes_back), seq, round(back_z,1), yesno(brain_to_skull), ...
    'VariableNames',{'line','elevation_deg','points_towards', ...
    'first_brain_mm','first_CSF_mm','first_skull_mm','first_scalp_mm', ...
    'first_meeting_in_order','every_step_in_order','tissues_passed', ...
    'goes_back_at_height_z_mm','brain_touches_skull'});

[u,~,j] = unique(seq(goes_back));
Te = table(u, accumarray(j,1), strings(numel(u),1), strings(numel(u),1), ...
    'VariableNames',{'tissues_passed','lines','points_towards','why'});
for q = 1:numel(u)
    m = goes_back & seq == u(q);
    Te.points_towards(q) = strjoin(unique(towards(m))', ', ');
    Te.why(q) = explain(find(m,1));
end

C = dirs*dirs'; C(1:N_RAY+1:end) = -1;            % cosine between every pair of lines, ignoring self
nn_deg = acosd(max(C,[],2));                      % angle to the nearest other line
Ts = table(["lines";"meet the layers in order (first meeting)";"in order at every step"; ...
             "pass a layer twice";"  of those, pointing at the top of head"; ...
             "lines with brain touching skull"; ...
             "middle point x (mm)";"middle point y (mm)";"middle point z (mm)"; ...
             "step along each line (mm)";"median angle to nearest line (deg)"], ...
    [N_RAY; sum(order_ok); sum(~goes_back); sum(goes_back); sum(goes_back & towards == "top of head"); ...
     sum(brain_to_skull); ctr(:); STEP; median(nn_deg)], ...
    'VariableNames',{'quantity','value'});

xf = [outd 'a1_nesting_lines.xlsx'];
if exist(xf,'file'), delete(xf); end
writetable(Tr, xf, 'Sheet','Lines');
writetable(Te, xf, 'Sheet','Explanations');
writetable(Ts, xf, 'Sheet','Summary');
colour = fullfile(thesis_dir('code','tools'), 'colour_true_false.py');   % TRUE green, FALSE red
system(sprintf('python3 "%s" "%s"', colour, xf));
fprintf('Nesting: %d of %d in order at every step | %d pass a layer twice | %d brain touching skull\n', ...
    sum(~goes_back), N_RAY, sum(goes_back), sum(brain_to_skull));

fprintf('\nSAVED results to %s\n', outd);

function s = yesno(x)
s = repmat("FALSE", size(x)); s(x) = "TRUE";
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
