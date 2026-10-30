%A1_METHODS_FIGURES  Draw the head model, the three electrode layouts and the set-up of two validation checks.

% Aim served: Aim 1 (Methods figures).
%
% Inputs : models/baseline/a1_baseline_*.mat, models/montages/a1_montages.mat,
%          results/aim1/a1_analytic_validation.mat
% Outputs: figures/methods/m1_healthy_head.png, m2_electrode_layouts.png,
%          m6_validation_setup.png
% Imogen Ryan - Honours thesis.

clear; close all force;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
out  = fullfile(root,'figures','methods'); if ~exist(out,'dir'), mkdir(out); end

% Tissue colours, in the order scalp, skull, CSF, brain
tis  = ["scalp","skull","CSF","brain"];
tisLab = ["Scalp","Skull","CSF","Brain"];                % names as shown on the figure
col  = [0.93 0.78 0.66; 0.95 0.95 0.88; 0.45 0.70 0.95; 0.85 0.55 0.60];
skin = [0.90 0.80 0.70];

L   = load(fullfile(root,'models','baseline','a1_baseline_m19_1020.mat'));
img = L.img; info = L.info; fm = img.fwd_model;
nd  = fm.nodes*1000;                                   % node positions (mm)
sig = [info.sigma.scalp info.sigma.skull info.sigma.csf info.sigma.brain];   % conductivity (S/m)
Mn  = load(fullfile(root,'models','montages','a1_montages.mat'),'M');   % electrode names

%% Figure 1: the healthy head
f = figure('Units','centimeters','Position',[1 1 17 8],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
tl = tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');

% (a) Scalp surface
ax = nexttile(tl);
head_surface(ax, fm, skin, 1);
view(ax,[-60 15]); camlight(ax,'headlight');
title(ax,'(a) Scalp surface (SAH262)','FontWeight','normal');

% (b) Cut down the middle, coloured by tissue
ax = nexttile(tl); hold(ax,'on');
slice_tissues(ax, img, info, 0, tis, col);
leg = strings(1,4);
for k = 1:4, leg(k) = sprintf('%s  %.2f S/m', tisLab(k), sig(k)); end
h = gobjects(1,4);
for k = 1:4, h(k) = patch(ax,NaN,NaN,col(k,:)); end
legend(ax, h, leg, 'Location','southeast','Box','off');
text(ax,-100,-60,'Front'); text(ax,85,-60,'Back');
axis(ax,'equal','off');
title(ax,'(b) Cut through the middle: four layers','FontWeight','normal');
exportgraphics(f, fullfile(out,'m1_healthy_head.png'),'Resolution',300);

%% Figure 2: the three electrode layouts, seen from above
lay  = ["m19_1020","m32_1010","m16_twoplane"];
ttl  = ["(a) 19 electrodes, 10-20","(b) 32 electrodes, 10-10","(c) 16 electrodes, two rings"];
f = figure('Units','centimeters','Position',[1 1 17 17],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
tl = tiledlayout(f,2,4,'Padding','compact','TileSpacing','compact');
tile = [1 3 6];                                        % (a) and (b) on the top row, (c) centred below
for m = 1:3
    B  = load(fullfile(root,'models','baseline',"a1_baseline_"+lay(m)+".mat"));
    ec = B.info.elec_centre_m*1000;                    % electrode centres (mm)
    ax = nexttile(tl, tile(m), [1 2]); hold(ax,'on');
    head_surface(ax, B.img.fwd_model, skin, 1);
    plot3(ax, ec(:,1), ec(:,2), ec(:,3)+2, 'o', 'MarkerSize',7, ...
          'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
    nm = string(Mn.M(strcmp({Mn.M.name},lay(m))).labels);
    ztop = max(B.img.fwd_model.nodes(:,3))*1000 + 20;     % labels above the scalp so the head never hides them
    text(ax, ec(:,1), ec(:,2)+8, ztop*ones(size(ec,1),1), nm, 'HorizontalAlignment','center');
    text(ax, 0, -118, ztop, 'nose', 'HorizontalAlignment','center');
    view(ax,[0 90]); camlight(ax,'headlight'); axis(ax,'equal','off');
    ndm = B.img.fwd_model.nodes*1000;
    xlim(ax, [min(ndm(:,1))-5 max(ndm(:,1))+5]); ylim(ax, [min(ndm(:,2))-5 max(ndm(:,2))+15]);   % head fills the tile, room for the title
    title(ax, sprintf('%s  (%d readings)', ttl(m), numel(fwd_solve(B.img).meas)), 'FontWeight','normal');
end
exportgraphics(f, fullfile(out,'m2_electrode_layouts.png'),'Resolution',300);

%% Figure 6: what the head looked like during two validation checks
Gv = load(fullfile(root,'results','aim1','a1_analytic_validation.mat'));
ctr_ball = Gv.CTR*1000;                                % test ball centre (mm)
f = figure('Units','centimeters','Position',[1 1 17 8],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
tl = tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');

% (a) Nesting check: straight lines out from the middle of the brain
ax = nexttile(tl); hold(ax,'on');
slice_tissues(ax, img, info, 0, tis, col);
c0 = mean(nd,1);                                       % middle of the head (mm)
th = linspace(0,2*pi,25); th(end) = [];
for k = 1:numel(th)
    d = [0 cos(th(k)) sin(th(k))];
    p = c0 + [0;131.8]*d;
    plot(ax, p(:,2), p(:,3), 'k-', 'LineWidth',0.7);
end
plot(ax, c0(2), c0(3), 'ko','MarkerFaceColor','y','MarkerSize',7);
axis(ax,'equal','off');
title(ax,'(a) Layer test: 400 lines out from the centre','FontWeight','normal');

% (b) Sphere check: a small ball in the brain, radius 3 to 20 mm
ax = nexttile(tl); hold(ax,'on');
slice_tissues(ax, img, info, ctr_ball(1), tis, col);
r = [3 4 5 6 8 10 13 16 20];                           % test ball radii (mm)
for k = 1:numel(r)
    plot(ax, ctr_ball(2)+r(k)*cos(linspace(0,2*pi,100)), ...
         ctr_ball(3)+r(k)*sin(linspace(0,2*pi,100)), 'k-', 'LineWidth',0.8);
end
axis(ax,'equal','off');
title(ax,'(b) Sphere test: a ball in the brain','FontWeight','normal');
exportgraphics(f, fullfile(out,'m6_validation_setup.png'),'Resolution',300);
fprintf('Saved Methods figures to %s\n', out);


function head_surface(ax, fm, c, a)
% Outer surface of the head, from the boundary faces
nd = fm.nodes*1000;
F  = fm.boundary;
% Turn every face to point outwards, so the lighting is even
n  = cross(nd(F(:,2),:)-nd(F(:,1),:), nd(F(:,3),:)-nd(F(:,1),:), 2);
fc = (nd(F(:,1),:)+nd(F(:,2),:)+nd(F(:,3),:))/3 - mean(nd,1);
in = sum(n.*fc,2) < 0;
F(in,:) = F(in,[1 3 2]);
patch(ax,'Faces',F,'Vertices',nd,'FaceColor',c,'EdgeColor','none','FaceAlpha',a, ...
      'FaceLighting','gouraud','AmbientStrength',0.5,'DiffuseStrength',0.6,'SpecularStrength',0);
axis(ax,'equal','off');
end

function slice_tissues(ax, img, info, x0, tis, col)
% Cut the head at x = x0 (mm) and colour each piece by its tissue
im = img; im.elem_data = zeros(numel(info.ref),1);
for k = 1:4, im.elem_data(info.ref == region_index(info.names,tis(k))) = k; end
s  = mdl_slice_mesher(im, [x0/1000 inf inf]);
p  = s.fwd_model.nodes*1000;
patch(ax,'Faces',s.fwd_model.elems,'Vertices',p(:,2:3),'FaceVertexCData',col(round(s.elem_data),:), ...
      'FaceColor','flat','EdgeColor','none');                % plotted as (front-back, up-down)
end
