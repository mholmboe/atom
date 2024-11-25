%% atom_rotate.m
% * This function rotate the atom randomly or by the angles given by the
% rotate vector
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = rotate_atom(atom,Box_dim)
% # atom = rotate_atom(atom,Box_dim,'random')
% # atom = rotate_atom(atom,Box_dim,[0 0 90])
% # atom = rotate_atom(atom,Box_dim,[0 0 90],ind)

function atom = rotate_atom(atom,Box_dim,varargin)

if nargin>2
    rotate=varargin{1};
else
    rotate='random';
end

if nargin>3
    ind=varargin{2};
else
    ind=0;
end

angles=[0 0 0];
if ischar(rotate)
    angles=[180-360*rand(1) 180-360*rand(1) 180-360*rand(1)];
elseif iscell(rotate)
    if size(rotate,2)==3
        for i=1:3
            if iscell(rotate(i))
                angles(i)=90*rand(i);
            else
                angles(i)=cell2mat(rotate(i));
            end
        end
    end
else
    if size(rotate,2)==3
%         disp('Set the angles to..')
%         rotate
        angles(1)=rotate(1);
        angles(2)=rotate(2);
        angles(3)=rotate(3);
    else
        disp('Did not catch any angles')
    end
end

alfa=angles(1);
beta=angles(2);
gamma=angles(3);

if ind==0
    atom = COM_atom(atom,Box_dim);
else
    COM=[atom(ind).x atom(ind).y atom(ind).z];
end

x_shift=num2cell([atom.x]-COM(1)); [atom.x]=deal(x_shift{:});
y_shift=num2cell([atom.y]-COM(2)); [atom.y]=deal(y_shift{:});
z_shift=num2cell([atom.z]-COM(3)); [atom.z]=deal(z_shift{:});

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
XYZ_data=XYZ_data*roty(alfa)*rotx(beta)*rotz(gamma);

x_rot=num2cell(XYZ_data(:,1)+COM(1)); [atom.x]=deal(x_rot{:});
y_rot=num2cell(XYZ_data(:,2)+COM(2)); [atom.y]=deal(y_rot{:});
z_rot=num2cell(XYZ_data(:,3)+COM(3)); [atom.z]=deal(z_rot{:});


try atom = rmfield(atom,'Mw'); catch; end
try atom = rmfield(atom,'element'); catch; end
try
    atom = rmfield(atom,'COM_x');
    atom = rmfield(atom,'COM_y');
    atom = rmfield(atom,'COM_z');
catch
end

% assignin('caller','rot_atom',rot_atom);

% write_atom_gro(atom,Box_dim,'out.gro');
%
% plot_atom('out.gro')

end

function rotmat = rotx(alpha)
% rotate in the direction of y->z, counter-clockwise
rotmat = [1 0 0;0 cosd(alpha) -sind(alpha); 0 sind(alpha) cosd(alpha)];
end

function rotmat = roty(beta)
% rotate in the direction of z->x, counter-clockwise
rotmat = [cosd(beta) 0 sind(beta); 0 1 0; -sind(beta) 0 cosd(beta)];
end

function rotmat = rotz(gamma)
% rotate in the direction of x->y, counter-clockwise
rotmat = [cosd(gamma) -sind(gamma) 0; sind(gamma) cosd(gamma) 0; 0 0 1];
end
