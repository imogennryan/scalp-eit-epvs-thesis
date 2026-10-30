%A2_METHODS_FIGURES  Draw the scalp points, the depths along one line and the effect of each fluid amount on conductivity.

% Aim served: Aim 2 (Methods figures).
%
% Inputs : models/baseline/a1_baseline_m19_1020.mat, models/montages/a1_montages.mat,
%          a2_place_perturbation, a2_mix_conductivity_L, a2_reachable_depth
% Outputs: figures/methods/m3_scalp_points.png, m4_fluid_amounts.png, m5_depths.png
% Imogen Ryan - Honours thesis.

clear; close all force;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'code','aim1_head_model'));
out  = fullfile(root,'figures','methods'); if ~exist(out,'dir'), mkdir(out); end

B    = load(fullfile(root,'models','baseline','a1_baseline_m19_1020.mat'));
img  = B.img; info = B.info; fm = img.fwd_model;
nd   = fm.nodes*1000;                                   % node positions (mm)
ctr  = mean(nd,1);                                      % head centre (mm)
Mn   = load(fullfile(root,'models','montages','a1_montages.mat'),'pos','names');

pts    = ["Cz","Fz","C3","Pz","T7","O1"];               % scalp points used in the sweep
depths = 25:5:60;                                       % target depths (mm)
burden = [0.0025 0.005 0.01 0.02 0.05 0.10];            % fluid share of the region
D_MM   = 20;                                            % region diameter (mm)
tis    = ["scalp","skull","CSF","brain"];
col    = [0.93 0.78 0.66; 0.95 0.95 0.88; 0.45 0.70 0.95; 0.85 0.55 0.60];
skin   = [0.90 0.80 0.70];

% Scalp positions and the deepest depth each line can reach
xyz = zeros(numel(pts),3); dmax = zeros(numel(pts),1);
for k = 1:numel(pts)
    xyz(k,:) = Mn.pos(strcmp(Mn.names,pts(k)),:);       % mm
    dmax(k)  = a2_reachable_depth(img, info, xyz(k,:)/1000);
end

%% Figure 3: the six scalp points and their lines into the head
f = figure('Units','centimeters','Position',[1 1 17 8],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
tl = tiledlayout(f,1,2,'Padding','compact','TileSpacing','compact');
vw = {[90 10], [0 90]};
ttl = ["(a) From the left side","(b) From above, nose at the bottom"];
for v = 1:2
    ax = nexttile(tl); hold(ax,'on');
    head_surface(ax, fm, skin, 0.35);
    for k = 1:numel(pts)
        u = (ctr - xyz(k,:)) / norm(ctr - xyz(k,:));    % direction into the head
        p = xyz(k,:) + [0; dmax(k)]*u;
        plot3(ax, p(:,1), p(:,2), p(:,3), 'k-', 'LineWidth',1.5);
        plot3(ax, xyz(k,1), xyz(k,2), xyz(k,3), 'o', 'MarkerSize',9, ...
              'MarkerFaceColor',[0.85 0.1 0.1], 'MarkerEdgeColor','k');
        text(ax, xyz(k,1)*1.08, xyz(k,2)*1.08, xyz(k,3)*1.08+4, pts(k), ...
             'FontWeight','bold','FontSize',10,'HorizontalAlignment','center');
    end
    view(ax, vw{v}); camlight(ax,'headlight'); axis(ax,'equal','off');
    if v == 2, ylim(ax,[-130 130]); end
    title(ax, ttl(v), 'FontWeight','normal');
end
exportgraphics(f, fullfile(out,'m3_scalp_points.png'),'Resolution',300);

%% Figure 5: the depths along the line below Cz
k  = 1;                                                 % Cz
u  = (ctr - xyz(k,:)) / norm(ctr - xyz(k,:));
f = figure('Units','centimeters','Position',[1 1 12 10],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
ax = axes(f); hold(ax,'on');
slice_tissues(ax, img, info, 0, tis, col);
p = xyz(k,:) + [0; dmax(k)]*u;
plot(ax, p(:,2), p(:,3), 'k--', 'LineWidth',1);
opt = struct('mode','radial','diam_mm',D_MM,'anchor_xyz_m',xyz(k,:)/1000);
cm  = flipud(parula(numel(depths)+2)); cm = cm(2:end-1,:);   % light = shallow
for j = 1:numel(depths)
    opt.depth_mm = depths(j);
    try
        [~, pl] = a2_place_perturbation(img, info, opt);
    catch
        continue                                        % region does not fit at this depth
    end
    if depths(j) + D_MM/2 > dmax(k), continue; end      % same fit rule as the sweep
    c = pl.centre_mm;
    th = linspace(0,2*pi,100);
    plot(ax, c(2)+D_MM/2*cos(th), c(3)+D_MM/2*sin(th), '-', 'Color',cm(j,:), 'LineWidth',1.8);
    if mod(j,2) == 1                                     % labels alternate right and left so they do not overlap
        text(ax, c(2)+13, c(3), sprintf('%d mm', depths(j)));
    else
        text(ax, c(2)-13, c(3), sprintf('%d mm', depths(j)), 'HorizontalAlignment','right');
    end
end
plot(ax, xyz(k,2), xyz(k,3), 'o','MarkerSize',9,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
text(ax, xyz(k,2), xyz(k,3)+8, 'Cz', 'FontWeight','bold','HorizontalAlignment','center');
axis(ax,'equal','off'); ylim(ax,[-60 100]);
exportgraphics(f, fullfile(out,'m5_depths.png'),'Resolution',300);

%% Figure 4: what each amount of fluid does to the conductivity, for round pockets and tubes
opt.depth_mm = 40;                                      % 40 mm below Cz
[frac, pl] = a2_place_perturbation(img, info, opt);
c  = pl.centre_mm;
sb = info.sigma.brain; si = info.sigma.csf;            % brain and CSF (S/m)
Ls   = [1/3 0];                                         % shape number of the pockets
Llab = {{'Round','pockets,','L = 1/3'}, {'Tubes along','the current,','L = 0'}};   % row labels, one tile each
f = figure('Units','centimeters','Position',[1 1 17 7.5],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width
tl = tiledlayout(f,2,numel(burden)+1,'Padding','compact','TileSpacing','compact');
for r = 1:numel(Ls)
  ax = nexttile(tl); axis(ax,'off'); xlim(ax,[0 1]); ylim(ax,[0 1]);
  text(ax, 0.5, 0.5, Llab{r}, 'HorizontalAlignment','center');
  for b = 1:numel(burden)
    fluid = full(frac)*burden(b);                       % fluid share of each element
    im = img; m = fluid > 0;
    im.elem_data(m) = a2_mix_conductivity_L(sb, si, fluid(m), Ls(r));
    s  = mdl_slice_mesher(im, [c(1)/1000 inf inf]);
    q  = s.fwd_model.nodes*1000;
    ax = nexttile(tl); hold(ax,'on');
    patch(ax,'Faces',s.fwd_model.elems,'Vertices',q(:,2:3),'FaceVertexCData',s.elem_data, ...
          'FaceColor','flat','EdgeColor',[0.6 0.6 0.6],'LineWidth',0.1);
    th = linspace(0,2*pi,100);
    plot(ax, c(2)+D_MM/2*cos(th), c(3)+D_MM/2*sin(th), 'k--');
    clim(ax,[sb 0.27]); colormap(ax, parula);
    axis(ax,'equal','off'); xlim(ax, c(2)+[-18 18]); ylim(ax, c(3)+[-18 18]);
    title(ax, {sprintf('%g%%', burden(b)*100), sprintf('%.1f mm^3', burden(b)*pi/6*D_MM^3)}, ...
          'FontWeight','normal');
  end
end
cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = 'Conductivity (S/m)';
exportgraphics(f, fullfile(out,'m4_fluid_amounts.png'),'Resolution',300);
fprintf('Saved Methods figures to %s\n', out);
fprintf('Deepest reachable depth (mm): %s\n', strjoin(compose('%s %.1f', pts', dmax), ', '));


function head_surface(ax, fm, c, a)
% Outer surface of the head, every face turned to point outwards
nd = fm.nodes*1000;
F  = fm.boundary;
n  = cross(nd(F(:,2),:)-nd(F(:,1),:), nd(F(:,3),:)-nd(F(:,1),:), 2);
fc = (nd(F(:,1),:)+nd(F(:,2),:)+nd(F(:,3),:))/3 - mean(nd,1);
in = sum(n.*fc,2) < 0;
F(in,:) = F(in,[1 3 2]);
patch(ax,'Faces',F,'Vertices',nd,'FaceColor',c,'EdgeColor','none','FaceAlpha',a, ...
      'FaceLighting','gouraud','AmbientStrength',0.5,'DiffuseStrength',0.6,'SpecularStrength',0);
end

function slice_tissues(ax, img, info, x0, tis, col)
% Cut the head at x = x0 (mm) and colour each piece by its tissue
im = img; im.elem_data = zeros(numel(info.ref),1);
for k = 1:4, im.elem_data(info.ref == region_index(info.names,tis(k))) = k; end
s  = mdl_slice_mesher(im, [x0/1000 inf inf]);
p  = s.fwd_model.nodes*1000;
patch(ax,'Faces',s.fwd_model.elems,'Vertices',p(:,2:3),'FaceVertexCData',col(round(s.elem_data),:), ...
      'FaceColor','flat','EdgeColor','none');
end
