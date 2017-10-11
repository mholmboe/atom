%% atom_rotate.m
% * This function rotate the atom
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se


%% Examples
% * atom = rotate_atom(atom,Box_dim,'random')
% * atom = rotate_atom(atom,Box_dim,[0 0 90])

function atom = rotate_atom(atom,Box_dim,rotate)


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
        disp('Set the angles to..')
        angles(1)=rotate(1);
        angles(2)=rotate(2);
        angles(3)=rotate(3);
    else
        disp('Did not catch any angles')
    end
    
end

alfa=angles(1)
beta=angles(2)
gamma=angles(3)

atom = COM_atom(atom,Box_dim);
x_shift=num2cell([atom.x]-COM(1)); [atom.x]=deal(x_shift{:});
y_shift=num2cell([atom.y]-COM(2)); [atom.y]=deal(y_shift{:});
z_shift=num2cell([atom.z]-COM(3)); [atom.z]=deal(z_shift{:});

XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
XYZ_data=XYZ_data*roty(alfa)*rotx(beta)*rotz(gamma);

x_rot=num2cell(XYZ_data(:,1)+COM(1)); [atom.x]=deal(x_rot{:});
y_rot=num2cell(XYZ_data(:,2)+COM(2)); [atom.y]=deal(y_rot{:});
z_rot=num2cell(XYZ_data(:,3)+COM(3)); [atom.z]=deal(z_rot{:});

rot_atom=atom;
atom = rmfield(atom,'Mw');
atom = rmfield(atom,'element');
atom = rmfield(atom,'COM_x');
atom = rmfield(atom,'COM_y');
atom = rmfield(atom,'COM_z');

assignin('caller','rot_atom',rot_atom);

% write_atom_gro(atom,Box_dim,'out.gro');
%
% plot_atom('out.gro')

