%% addvsites_atom.m
% * This function generates virtual (ghost ) sites on top of real particle
% sites. Make sure to pass on the full atom struct in the function call.
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # new_atom=addvsites_atom(atom,'Na','NaV') % Basic input arguments, Real atomtype, virtual atomtype, virtual residue name
% # new_atom=addvsites_atom(atom,'Na','NaV',[0 0 10]) % Translates the copied sites
% # new_atom=addvsites_atom(atom,'Na','NaV'[0 0 10],10) % Translates and copies only 10 sites
%
function atom = addvsites_atom(atom,atomtype,new_atomtype,varargin)

% Find all original atomtypes to replace
ind_atomtype=find(strcmpi([atom.type],atomtype));

if nargin>4
    trans_vec=varargin{1};
else
    trans_vec=[0 0 0];
end

% Duplicate the atomtype entries into a new In_atom struct
if nargin == 5
    % Uses the num last entries...
    num=cell2mat(varargin(2));
    virtual_atom=atom(ind_atomtype(end-num+1:end));
else
    % Uses all entries...
    virtual_atom=atom(ind_atomtype);
end

% Rename all
[virtual_atom.type]=deal({new_atomtype});
[virtual_atom.fftype]=deal({new_atomtype});

if max(abs(trans_vec))>0
    virtual_atom = translate_atom(virtual_atom,trans_vec,'all')
end

for i=1:numel(ind_atomtype)
    [virtual_atom(i).index]=[virtual_atom(i).index]+0.5;
end

atom=[atom virtual_atom];

[ind_val,ind]=sort([atom.index]);

atom=update_atom(atom(ind));

end