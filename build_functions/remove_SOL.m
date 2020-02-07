%% remove_SOL.m
% * This function removes water molecules between the lo and hi limits
% * along the dim dimension. Does it work for tip4p and tip5p water?
% * Untested for quite some time
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = remove_SOL(atom,atomname,0,10,'z')
% # atom = remove_SOL(atom,atomname,0,10,'z','OW')

function atom = remove_SOL(atom,atomname,lo,hi,dim,varargin)
%%

nAtoms=size(atom,2);

ind=find(strncmpi([atom.type],atomname,2));

if strcmp(dim,'x')
    atomcoords=[atom.x];
elseif strcmp(dim,'y')
    atomcoords=[atom.y];
elseif strcmp(dim,'z')
    atomcoords=[atom.z];
end

ind_lo=find(atomcoords>lo);ind_lo=intersect(ind_lo,ind);
ind_hi=find(atomcoords<hi);ind_hi=intersect(ind_hi,ind);
ind=intersect(ind_lo,ind_hi);
if strncmpi(atomname,'OW',2)
    ind_HW1=ind+1;ind_HW2=ind+2;
    ind=sort([ind ind_HW1 ind_HW2]);
end

if nargin == 6
    if strncmpi(atomname,'OW',2)
        nrmSOL=3*cell2mat(varargin(1));
    else
        nrmSOL=cell2mat(varargin(1));
    end
else
    nrmSOL=length(ind);
end

length(ind)

ind=ind(1:nrmSOL);

change_ind=sort(setdiff(1:nAtoms,ind));
atom=atom(change_ind);
atom=update_atom(atom);
% %