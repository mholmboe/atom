%% sort_atom.m
% * This function reorders the atoms in an atom struct. Useful for instance
% * for creating united-atom structures from all-atom structures, by simply
% * ignoring the non-polar H's from the original list of atoms, or
% * reordering the atom struct with respect to residue name or atom type
% *
% * In case of reordering the the atom order, neworder is a [1xn] array of
% * n index numbers with a new order. Else if varargin is 'resname' or
% * 'atomtype', neworder is a cell list of 'stringnames'
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = reorder_atom(atom,data_val) % Orders according to index
%

function atom_sorted = sort_atom(atom,data_val,varargin)

% Determine sort direction
orderDir = 'ascend';
if nargin>2 && ischar(varargin{1}) && strcmpi(varargin{1}, 'descend')
    orderDir = 'descend';
end


% Total number of atoms
nAtoms = size(atom,2);

% Number of groups
Groups = unique([atom.type],'stable');
nGroups= numel(Groups);

% Identify unique atom types in stable order
types = [atom.type];
Groups = unique(types, 'stable');
nGroups = numel(Groups);

% Preallocate sorted struct
atom_sorted = atom;

% Loop over each type group and sort by data_val
for g = 1:nGroups
    % find indices for this type
    idx = find(strcmp(types, Groups{g}));
    % extract values for sorting
    vals = data_val(idx);
    % get sorting order within the group
    [~, sortIdx] = sort(vals, orderDir);
    % reorder atoms in this group
    atom_sorted(idx) = atom(idx(sortIdx));
end

end


