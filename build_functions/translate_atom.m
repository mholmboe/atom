%% translate_atom.m
% * This function translates the resname by a vector
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = translate_molid(atom,[x y z])
% # atom = translate_molid(atom,[x y z],'all')
% # atom = translate_molid(atom,[x y z],'SOL')

function atom = translate_atom(atom,trans_vec,varargin)

nAtoms=size([atom.x],2);

if nargin>2
    Resname=varargin(1);
else
    Resname='all';
end

if strcmp(Resname,'all')
    ind_resname=1:nAtoms;
elseif strcmp(Resname,'All')
    ind_resname=1:nAtoms;
elseif strcmp(Resname,'ALL')
    ind_resname=1:nAtoms;
else
    ind_resname=find(strcmpi([atom.resname],Resname));
end

if size(trans_vec,1)==1
    x_shift=num2cell([atom(ind_resname).x]+trans_vec(1)); [atom(ind_resname).x]=deal(x_shift{:});
    
    y_shift=num2cell([atom(ind_resname).y]+trans_vec(2)); [atom(ind_resname).y]=deal(y_shift{:});
    
    z_shift=num2cell([atom(ind_resname).z]+trans_vec(3)); [atom(ind_resname).z]=deal(z_shift{:});
else
    x_shift=num2cell([atom(ind_resname).x]+trans_vec(:,1)');[atom(ind_resname).x]=deal(x_shift{:});
    
    y_shift=num2cell([atom(ind_resname).y]+trans_vec(:,2)'); [atom(ind_resname).y]=deal(y_shift{:});
    
    z_shift=num2cell([atom(ind_resname).z]+trans_vec(:,3)'); [atom(ind_resname).z]=deal(z_shift{:});    
end

% assignin('caller','atom',atom);
assignin('caller','XYZ_data',[[atom.x]' [atom.y]' [atom.z]']);
assignin('caller','XYZ_labels',[atom.type]');
