%% atom_spiral.m
% * This function spiral the atom randomly or by the angles given by the
% spiral vector
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = spiral_atom(atom,Box_dim,[0 0 90])
%
function atom = spiral_atom(atom,Box_dim,angles,varargin)

if nargin>3
    translate_vec=varargin{1};
else
    translate_vec=0;
end

if ~size(angles,2)==3
    disp('Did not catch any angles')
end

alfa=angles(1);
beta=angles(2);
gamma=angles(3);

if sum(abs(translate_vec))>0
    atom=translate_atom(atom,translate_vec);
end

x_shift=num2cell([atom.x]); [atom.x]=deal(x_shift{:});
y_shift=num2cell([atom.y]); [atom.y]=deal(y_shift{:});
z_shift=num2cell([atom.z]); [atom.z]=deal(z_shift{:});

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
XYZ_data=XYZ_data*roty(alfa)*rotx(beta)*rotz(gamma);

x_rot=num2cell(XYZ_data(:,1)); [atom.x]=deal(x_rot{:});
y_rot=num2cell(XYZ_data(:,2)); [atom.y]=deal(y_rot{:});
z_rot=num2cell(XYZ_data(:,3)); [atom.z]=deal(z_rot{:});

end

function rotmat = rotx(alpha)
% eml_assert_no_varsize(1,alpha);
% sigdatatypes.validateAngle(alpha,'rotx','ALPHA',{'scalar'});
% spiral in the direction of y->z, counter-clockwise
rotmat = [1 0 0;0 cosd(alpha) -sind(alpha); 0 sind(alpha) cosd(alpha)];
end

function rotmat = roty(beta)
% eml_assert_no_varsize(1,beta);
% sigdatatypes.validateAngle(beta,'roty','BETA',{'scalar'});
% spiral in the direction of z->x, counter-clockwise
rotmat = [cosd(beta) 0 sind(beta); 0 1 0; -sind(beta) 0 cosd(beta)];
end

function rotmat = rotz(gamma)
% eml_assert_no_varsize(1,gamma);
% sigdatatypes.validateAngle(gamma,'rotz','GAMMA',{'scalar'});
% spiral in the direction of x->y, counter-clockwise
rotmat = [cosd(gamma) -sind(gamma) 0; sind(gamma) cosd(gamma) 0; 0 0 1];
end
