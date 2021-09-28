%% remove_residues.m
% * This section is used to remove residues in the simulation box between
% * limits lo and hi in the dim dimension
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = remove_residues(atom,'SOL',10,20,'z') % Basic input arguments
%
function atom = remove_residues(atom,resnames,lo,hi,dim)

nAtoms=size(atom,2);

ind=zeros(1,size(atom,2));
for i = 1:length(resnames)
    ind=[ind find(strcmp([atom.resname],resnames(i)))];
end
ind=sort(ind);

atom = median_atom_func(atom);

if strcmpi(dim,'x')
    atomcoords=[atom.med_x];
elseif strcmpi(dim,'y')
    atomcoords=[atom.med_y];
elseif strcmpi(dim,'z')
    atomcoords=[atom.med_z];
end

ind_lo=find(atomcoords>lo);ind_lo=intersect(ind_lo,ind)
ind_hi=find(atomcoords<hi);ind_hi=intersect(ind_hi,ind)
ind=intersect(ind_lo,ind_hi);
if strcmpi(resname,'OW')
    ind_HW1=ind+1;ind_HW2=ind+2;
    ind=sort([ind ind_HW1 ind_HW2]);
end
change_ind=sort(setdiff(1:nAtoms,ind));
atom=atom(change_ind);
atom=atom_update(atom);
% %

