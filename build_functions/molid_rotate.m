%% molid_rotate.m
% * This function rotates the atom molid randomly
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = molid_rotate(atom,Box_dim,[1 2 4],'xy')
%
function atom = molid_rotate(atom,Box_dim,MolID,rotate_dim)

x_vec=[1 0 0];
y_vec=[0 1 0];
z_vec=[0 0 1];

for i=MolID
    
    molid_ind=ismember([atom.molid],i);
    rot_atom = atom(molid_ind);
    rot_atom = unwrap_atom_func(rot_atom,Box_dim,'xyz');
    
    princaxisangle=[max([rot_atom.x])-min([rot_atom.x]) max([rot_atom.y])-min([rot_atom.y]) max([rot_atom.z])-min([rot_atom.z])];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strfind(rotate_dim,'x');
        angle_x=rad2deg(atan2(norm(cross(princaxisangle,x_vec)),dot(princaxisangle,x_vec)))
    else
        angle_x=0;
    end
    if strfind(rotate_dim,'y');
        angle_y=rad2deg(atan2(norm(cross(princaxisangle,y_vec)),dot(princaxisangle,y_vec)))
    else
        angle_y=0;
    end
    if strfind(rotate_dim,'z');
        angle_z=rad2deg(atan2(norm(cross(princaxisangle,z_vec)),dot(princaxisangle,z_vec)))
    else
        angle_z=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strfind(rotate_dim,'-x');
        angle_x=rad2deg(atan2(norm(cross(princaxisangle,-x_vec)),dot(princaxisangle,x_vec)))
    end
    if strfind(rotate_dim,'-y');
        angle_y=rad2deg(atan2(norm(cross(princaxisangle,-y_vec)),dot(princaxisangle,y_vec)))
    end
    if strfind(rotate_dim,'-z');
        angle_z=rad2deg(atan2(norm(cross(princaxisangle,-z_vec)),dot(princaxisangle,z_vec)))
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strncmpi('random',rotate_dim,4);
        disp('Random rotation')
        angle_x=rand(1)*360;
        angle_y=rand(1)*360;
        angle_z=rand(1)*360;
    end
    
    
    rot_atom = COM_atom_func(rot_atom,Box_dim,i);
    x_shift=num2cell([rot_atom.x]-COM_molid(1)); [rot_atom.x]=deal(x_shift{:});
    y_shift=num2cell([rot_atom.y]-COM_molid(2)); [rot_atom.y]=deal(y_shift{:});
    z_shift=num2cell([rot_atom.z]-COM_molid(3)); [rot_atom.z]=deal(z_shift{:});
    
    XYZ_data=[[rot_atom.x]' [rot_atom.y]' [rot_atom.z]'];
    XYZ_data=XYZ_data*rotx(-angle_x)*roty(-angle_y)*rotz(-angle_z);
    
    x_rot=num2cell(XYZ_data(:,1)+COM_molid(1)); [atom(molid_ind).x]=deal(x_rot{:});
    y_rot=num2cell(XYZ_data(:,2)+COM_molid(2)); [atom(molid_ind).y]=deal(y_rot{:});
    z_rot=num2cell(XYZ_data(:,3)+COM_molid(3)); [atom(molid_ind).z]=deal(z_rot{:});
    
end

end

function rotmat = rotx(alpha)
eml_assert_no_varsize(1,alpha);
sigdatatypes.validateAngle(alpha,'rotx','ALPHA',{'scalar'});
% rotate in the direction of y->z, counter-clockwise
rotmat = [1 0 0;0 cosd(alpha) -sind(alpha); 0 sind(alpha) cosd(alpha)];
end

function rotmat = roty(beta)
eml_assert_no_varsize(1,beta);
sigdatatypes.validateAngle(beta,'roty','BETA',{'scalar'});
% rotate in the direction of z->x, counter-clockwise
rotmat = [cosd(beta) 0 sind(beta); 0 1 0; -sind(beta) 0 cosd(beta)];
end

function rotmat = rotz(gamma)
eml_assert_no_varsize(1,gamma);
sigdatatypes.validateAngle(gamma,'rotz','GAMMA',{'scalar'});
% rotate in the direction of x->y, counter-clockwise
rotmat = [cosd(gamma) -sind(gamma) 0; sind(gamma) cosd(gamma) 0; 0 0 1];
end


