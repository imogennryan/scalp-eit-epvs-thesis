function idx = region_index(names, want)

%REGION_INDEX  Find the materials.ref index of a tissue by its name.

% Aim served: Aim 1 and Aim 2 (shared helper).
%
% The order of the tissues in materials.names is not the same on every
% Tizzard mesh, so a tissue must be found by name, never by a fixed number.
%
%   idx = region_index(names, want)
%     names  region names as a cell array, e.g. strtrim(cellstr(S.materials.names))
%     want   tissue to find, 'scalp', 'skull', 'CSF' or 'brain' (not case sensitive)
%     idx    index into names, to compare against materials.ref
% Gives an error if the name is missing or appears more than once.
% Imogen Ryan - Honours thesis.

idx = find(strcmpi(strtrim(names), want));
if numel(idx) ~= 1
    error('region_index:notFound', ...
        'Region "%s" not found (or ambiguous) in [%s]', want, strjoin(names(:)',', '));
end
end
