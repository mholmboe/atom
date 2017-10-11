%% center_atom.m
% * This function centers the atom with respect to the resname molecule. newBox_dim could be a new or old Box_dim
% * atom is the atom struct
% * Box_dim is the box dimension vector
% * resname is the Resnames ofthe atoms ou want to move, can be 'all'
% * dim is a string containing for example x, xy, xyz
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom = center_atom(atom,Box_dim,'all','xy')
% * atom = center_atom(atom,Box_dim,'SOL','Z')


function atom = center_atom(atom,Box_dim,varargin)
disp('Center to box')

nAtoms=size(atom,2);

if nargin == 2
    resname='all';
    dim='xyz';
elseif nargin == 3
    resname=varargin{1};
    dim='xyz';
elseif nargin == 4
    resname=varargin{1};
    dim=varargin{2};
end


% Box_dim=newBox_dim;

% atom = wrap_atom(atom,Box_dim);

if strcmpi(resname,'all')
    ind_resname=1:nAtoms;
else
    if sum(find(strcmp([atom.resname],resname)))>0
        ind_resname=find(strcmp([atom.resname],resname));
    elseif sum(find(strcmp([atom.type],resname)))>0
        ind_resname=find(strcmp([atom.type],resname));
    end
end

if strfind(dim,'x') | strfind(dim,'X')
    disp('centering along x')
    x_shift=num2cell([atom.x]-median([atom(ind_resname).x])+Box_dim(1)/2); [atom(:).x]=deal(x_shift{:});
    x_shift=num2cell([atom.x]-mean([atom(ind_resname).x])+Box_dim(1)/2); [atom(:).x]=deal(x_shift{:});
end
if strfind(dim,'y') | strfind(dim,'Y')
    disp('centering along y')
    y_shift=num2cell([atom.y]-median([atom(ind_resname).y])+Box_dim(2)/2); [atom(:).y]=deal(y_shift{:});
    y_shift=num2cell([atom.y]-mean([atom(ind_resname).y])+Box_dim(2)/2); [atom(:).y]=deal(y_shift{:});
end
if strfind(dim,'z') | strfind(dim,'Z')
    disp('centering along z')
    median([atom(ind_resname).z])
    z_shift=num2cell([atom.z]-median([atom(ind_resname).z])+Box_dim(3)/2);   [atom(:).z]=deal(z_shift{:});
    z_shift=num2cell([atom.z]-mean([atom(ind_resname).z])+Box_dim(3)/2); [atom(:).z]=deal(z_shift{:});
end

assignin('caller','Box_dim',Box_dim);
assignin('caller','XYZ_data',[[atom.x]' [atom.y]' [atom.z]']);
assignin('caller','XYZ_labels',[atom.type]');
