%A2_METHODS_FIGURE_CSF  Draw a few 20 mm fluid regions in a cut through the head.

% Aim served: Aim 2 (Methods figure).
%
% Inputs : models/baseline/a1_baseline_m19_1020.mat, models/montages/a1_montages.mat,
%          a2_place_perturbation
% Outputs: figures/methods/m8_fluid_vs_csf.png
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

D_MM   = 20;                                            % region diameter (mm)
tis    = ["scalp","skull","CSF","brain"];
tisLab = ["Scalp","Skull","CSF","Brain"];                % names as shown on the figure
col    = [0.93 0.78 0.66; 0.95 0.95 0.88; 0.45 0.70 0.95; 0.85 0.55 0.60];

% A few configurations on the midline, so a cut down the middle shows them exactly
cfg = {"Fz",30; "Cz",25; "Cz",50; "Pz",45};             % scalp point, depth (mm)

f = figure('Units','centimeters','Position',[1 1 14 9],'Color','w','Visible','off'); theme(f,'light');
set(f, 'DefaultAxesFontName','Times New Roman', 'DefaultTextFontName','Times New Roman', 'DefaultLegendFontName','Times New Roman', 'DefaultColorbarFontName','Times New Roman', 'DefaultAxesFontSize',10, 'DefaultTextFontSize',10, 'DefaultLegendFontSize',10, 'DefaultColorbarFontSize',10);   % thesis font, 10 pt at the printed width

%% The regions drawn in the head
ax = axes(f); hold(ax,'on');
slice_tissues(ax, img, info, 0, tis, col);
th = linspace(0,2*pi,100);
for j = 1:size(cfg,1)
    xyz = Mn.pos(strcmp(Mn.names,cfg{j,1}),:);          % scalp point (mm)
    opt = struct('mode','radial','diam_mm',D_MM,'anchor_xyz_m',xyz/1000,'depth_mm',cfg{j,2});
    [~, pl] = a2_place_perturbation(img, info, opt);
    c = pl.centre_mm;
    u = (ctr - xyz)/norm(ctr - xyz);                    % direction into the head
    p = xyz + [0; cfg{j,2}]*u;
    plot(ax, p(:,2), p(:,3), 'k:', 'LineWidth',1);
    fill(ax, c(2)+D_MM/2*cos(th), c(3)+D_MM/2*sin(th), [0.1 0.1 0.6], ...
         'FaceAlpha',0.35, 'EdgeColor',[0.1 0.1 0.6], 'LineWidth',1.5);
    plot(ax, xyz(2), xyz(3), 'o','MarkerSize',8,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
    text(ax, c(2)+13, c(3), sprintf('%s, %d mm', cfg{j,1}, cfg{j,2}));
end
pts = unique(string(cfg(:,1)),'stable');
for k = 1:numel(pts)
    xyz = Mn.pos(strcmp(Mn.names,pts(k)),:);
    text(ax, xyz(2), xyz(3)+8, pts(k), 'FontWeight','bold','HorizontalAlignment','center');
end
for k = 1:4                                             % tissue key
    patch(ax, 100+[0 8 8 0], 100-11*k+[0 0 8 8], col(k,:), 'EdgeColor',[0.5 0.5 0.5]);
    text(ax, 111, 104-11*k, tisLab(k));
end
axis(ax,'equal','off'); xlim(ax,[-125 145]); ylim(ax,[-70 100]);

exportgraphics(f, fullfile(out,'m8_fluid_vs_csf.png'),'Resolution',300);
fprintf('Saved %s\n', fullfile(out,'m8_fluid_vs_csf.png'));


function slice_tissues(ax, img, info, x0, tis, col)
% Cut the head at x = x0 (mm) and colour each piece by its tissue
im = img; im.elem_data = zeros(numel(info.ref),1);
for k = 1:4, im.elem_data(info.ref == region_index(info.names,tis(k))) = k; end
s  = mdl_slice_mesher(im, [x0/1000 inf inf]);
p  = s.fwd_model.nodes*1000;
patch(ax,'Faces',s.fwd_model.elems,'Vertices',p(:,2:3),'FaceVertexCData',col(round(s.elem_data),:), ...
      'FaceColor','flat','EdgeColor','none');
end
