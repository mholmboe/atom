%% scale_atom.m
% * This function scales the coordinates in the atom struct
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = scale_atom(atom,Box_dim,[1 1 1.5],'all')
% # atom = scale_atom(atom,Box_dim,0.9 0.9 0.9,'SOL')
%
function atom = scale_atom(atom,Box_dim,scale_vec,varargin)

if nargin>3
    Resname=varargin{1};
else
    Resname='ALL';
end

disp('Scaling the coordinates')

nAtoms=size([atom.x],2);

if strncmpi(Resname,'ALL',3)
    ind_resname=1:nAtoms;
elseif strncmpi(Resname,'SYSTEM',3)
    ind_resname=1:nAtoms;
else
    ind_resname=find(strcmp([atom.resname],Resname));
end

x_shift=num2cell([atom(ind_resname).x]*scale_vec(1)); [atom(ind_resname).x]=deal(x_shift{:});

y_shift=num2cell([atom(ind_resname).y]*scale_vec(2)); [atom(ind_resname).y]=deal(y_shift{:});

z_shift=num2cell([atom(ind_resname).z]*scale_vec(3)); [atom(ind_resname).z]=deal(z_shift{:});

Box_dim(1:3)=Box_dim(1:3).*scale_vec;
if numel(Box_dim(1,1:end))==9
    Box_dim(6)=Box_dim(6)*scale_vec(1); % Tested by scaling - replicating triclinic box
    Box_dim(8)=Box_dim(8)*scale_vec(3); %
    Box_dim(9)=Box_dim(9)*scale_vec(2); %
end

atom=update_atom(atom);

assignin('caller','Box_dim',Box_dim);
