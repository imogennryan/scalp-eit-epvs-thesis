%A1_RECIPROCITY_TABLE  Save the reciprocity check on the 19-electrode head as tables.

% Aim served: Aim 1.
%
% Every adjacent pair drives in turn and every pair is read, giving a 19 x 19
% grid. Driving pair i and reading pair j must equal driving j and reading i.
%
% Inputs : a1_build_head_model, models/montages/a1_montages.mat
% Outputs: results/aim1/a1_reciprocity_grid.csv (19 x 19 readings, mV),
%          results/aim1/a1_reciprocity_examples.csv (worked examples)
% Imogen Ryan - Honours thesis.

clear; close all force;
addpath(fileparts(mfilename('fullpath')));
outd = thesis_dir('results','aim1');

% Build the head and solve with every pair driving and every pair read
[img, info] = a1_build_head_model();
ne = numel(img.fwd_model.electrode);
[stim, msel] = mk_stim_patterns(ne, 1, '{ad}', '{ad}', {'meas_current'}, info.opt.amp);
img.fwd_model.stimulation = stim; img.fwd_model.meas_select = msel;
v = fwd_solve(img);
M = reshape(v.meas, ne, ne);          % M(j,i): read pair j while driving pair i (V)

% Names of the pairs
Mt  = load([thesis_dir('models','montages') 'a1_montages.mat'], 'M');
lab = Mt.M(strcmp({Mt.M.name}, info.opt.montage)).labels;
pair = arrayfun(@(p) string(sprintf('%s-%s', lab{p}, lab{mod(p,ne)+1})), (1:ne)');

% The grid, in mV
G = array2table(M*1e3, 'VariableNames', cellstr(pair), 'RowNames', cellstr(pair));
writetable(G, [outd 'a1_reciprocity_grid.csv'], 'WriteRowNames', true);

% Worked examples: drive i, read j, then swapped
ex = [1 10; 1 15; 6 13];               % Fp1-Fp2/Cz-C4, Fp1-Fp2/Pz-P4, F4-F8/P7-P3
E = table();
for q = 1:size(ex,1)
    i = ex(q,1); j = ex(q,2);
    E = [E; table(pair(i), pair(j), M(j,i)*1e3, M(i,j)*1e3, abs(M(j,i)-M(i,j))*1e3, ...
        'VariableNames', {'current_through','meter_on','reading_mV', ...
        'reading_swapped_mV','difference_mV'})]; %#ok<AGROW>
end
writetable(E, [outd 'a1_reciprocity_examples.csv']);

mismatch = max(max(abs(M - M.'))) / max(abs(M(:)));
fprintf('%d pairs, %d readings, each compared with its swapped partner\n', ne, ne^2);
fprintf('largest mismatch: %.2e of the largest reading\n', mismatch);
disp(E)

function d = thesis_dir(varargin)
% Folder inside this repository
r = fileparts(mfilename('fullpath'));
while ~exist(fullfile(r,'README.md'),'file'), r = fileparts(r); end
d = [fullfile(r, varargin{:}) filesep];
end
