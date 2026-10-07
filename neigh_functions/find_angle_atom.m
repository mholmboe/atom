%% find_angle_atom.m
% * This function finds the rows of an Angle_index (or Bond_index) that involve
% * certain atomtypes, and returns their row indices. It is meant for picking out
% * for instance all Si-Ob-Si angles, or every angle involving both Ob and Si, from
% * the full Angle_index produced by bond_atom() or bond_angle_dihedral_atom().
% *
% * Pass the index columns only, i.e. Angle_index (extra columns are ignored, the
% * atom indices are taken from columns 1-3) or Bond_index(:,1:2) for bonds. A
% * matrix with two columns is treated as pairs, one with three or more as triplets.
% *
% * The matching mode decides what 'involving' means:
% *   'contains'  (default) the triplet holds at least one atom of each listed type.
% *               Repeat a type to ask for two of them, so {'Ob','Ob'} requires two.
% *   'only'      every atom of the triplet is one of the listed types
% *   'exact'     the types match the given pattern, e.g. {'Ob','Si','Op'}, with the
% *               two ends allowed to swap. Use '*' as a wildcard for any type.
% *   'center'    only the central atom, i.e. column 2, has to match
% *
% * The returned indices plug straight into show_atom_springs, see the examples.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # ind = find_angle_atom(atom,Angle_index,'Ob') % Every angle involving an Ob
% # ind = find_angle_atom(atom,Angle_index,{'Ob','Si'}) % Involving both Ob and Si
% # ind = find_angle_atom(atom,Angle_index,{'Ob','Si','Op'}) % Involving all three
% # ind = find_angle_atom(atom,Angle_index,{'Ob','Si'},'mode','only') % Only Ob and Si
% # ind = find_angle_atom(atom,Angle_index,{'Ob','Si','Ob'},'mode','exact') % Ob-Si-Ob
% # ind = find_angle_atom(atom,Angle_index,{'*','Ob','*'},'mode','exact') % Centred on Ob
% # ind = find_angle_atom(atom,Angle_index,{'Si'},'mode','center') % The same thing
% # ind = find_angle_atom(atom,Bond_index(:,1:2),{'Si','Ob'},'mode','only') % Si-Ob bonds
% # show_atom_springs(atom,Box_dim,'spring_index',Angle_index(ind,1:3))
%

function [ind,types] = find_angle_atom(atom,index_matrix,sel_types,varargin)
%%

% --- options (name/value) ---
p = struct('mode','contains');
for k=1:2:numel(varargin)
    fn = fieldnames(p);
    hit = find(strcmpi(fn,varargin{k}),1);
    if isempty(hit)
        error('find_angle_atom:UnknownOption','Unknown option "%s".', num2str(varargin{k}));
    end
    p.(fn{hit}) = varargin{k+1};
end

if ischar(sel_types) || isstring(sel_types), sel_types = {char(sel_types)}; end
sel_types = cellfun(@char, sel_types, 'UniformOutput', false);

ind=[]; types={};
if isempty(index_matrix), return; end

nc = size(index_matrix,2);
if nc < 2
    error('find_angle_atom:BadIndex','The index matrix needs at least two columns.');
end
if nc == 2
    cols = 1:2;                       % a pair list
elseif nc == 3 && any(mod(index_matrix(:,3),1) ~= 0)
    % A three-column matrix whose last column is not integer is a Bond_index,
    % where column 3 holds the bond distance rather than a third atom.
    cols = 1:2;
    disp('Third column is not integer, treating this as a Bond_index of pairs')
else
    cols = 1:3;                       % a triplet list, e.g. Angle_index
end

T = [atom.type];
idx = index_matrix(:,cols);
if any(idx(:) < 1) || any(idx(:) > numel(atom))
    error('find_angle_atom:OutOfRange','The index matrix points outside the atom struct.');
end
tt = reshape(T(idx), size(idx));       % types of each member, same shape as idx

keep = false(size(idx,1),1);
switch lower(p.mode)
    case 'contains'
        % Each listed type must be present, counting repeats.
        ureq = unique(sel_types);
        need = cellfun(@(u) sum(strcmp(sel_types,u)), ureq);
        for r = 1:size(idx,1)
            ok = true;
            for q = 1:numel(ureq)
                if sum(strcmp(tt(r,:), ureq{q})) < need(q), ok = false; break; end
            end
            keep(r) = ok;
        end
    case 'only'
        for r = 1:size(idx,1)
            keep(r) = all(ismember(tt(r,:), sel_types));
        end
    case 'exact'
        if numel(sel_types) ~= numel(cols)
            error('find_angle_atom:BadPattern', ...
                  'Mode ''exact'' needs a pattern with %d types.', numel(cols));
        end
        pat = sel_types;
        rev = pat(end:-1:1);
        for r = 1:size(idx,1)
            keep(r) = match_(tt(r,:),pat) || match_(tt(r,:),rev);
        end
    case 'center'
        if numel(cols) < 3
            error('find_angle_atom:NoCenter','Mode ''center'' needs triplets, not pairs.');
        end
        keep = ismember(tt(:,2), sel_types);
    otherwise
        error('find_angle_atom:BadMode', ...
              'Unknown mode "%s", use contains | only | exact | center.', p.mode);
end

ind   = find(keep);
types = tt(keep,:);

end  % main function


% ===== local helpers =====
function tf = match_(t, pat)
% Element-wise type match, where '*' in the pattern matches anything.
tf = true;
for k = 1:numel(pat)
    if strcmp(pat{k},'*'), continue; end
    if ~strcmp(t{k}, pat{k}), tf = false; return; end
end
end
