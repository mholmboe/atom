%% write_atom_cif.m
% * This function writes a basic 'P 1' cif file with fractional coordinates
% from the atom struct
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_cif(atom,Box_dim,filename_out)
%
function write_atom_cif(atom,Box_dim,filename_out)

% atom=wrap_atom(atom,Box_dim); % Do we need this

if regexp(filename_out,'.cif') ~= false
    filename_out = filename_out;
else
    filename_out = strcat(filename_out,'.cif');
end

if numel(Box_dim)<9
    if numel(Box_dim)==1
        Box_dim(1)=Box_dim(1);
        Box_dim(2)=Box_dim(1);
        Box_dim(3)=Box_dim(1);
    end
    Box_dim(4:9)=0;
end

disp('Assuming P1 space group. Box can still be triclinic')

Box_dim(Box_dim<0.00001&Box_dim>-0.00001)=0;

lx=Box_dim(1);
ly=Box_dim(2);
lz=Box_dim(3);
xy=Box_dim(6);
xz=Box_dim(8);
yz=Box_dim(9);

a=lx;
b=(ly^2+xy^2)^.5;
c=(lz^2+xz^2+yz^2)^.5;

if numel(find(Box_dim(4:9)))>0
    alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)));
    beta=rad2deg(acos(xz/c));
    gamma=rad2deg(acos(xy/b));
else
    alfa=90;
    beta=90;
    gamma=90;
end

%     a=Box_dim(1);
%     b=Box_dim(2);
%     c=Box_dim(3);
%     xy=Box_dim(6);
%     xz=Box_dim(8);
%     yz=Box_dim(9);
%     lx = a;
%     ly = (b^2-xy^2)^.5;
%     lz = (c^2 - xz^2 - yz^2)^0.5;
%     alfa=rad2deg(acos((ly*yz+xy*xz)/(b*c)))
%     beta=rad2deg(acos(xz/c));
%     gamma=rad2deg(acos(xy/b));


for i=1:size(atom,2)
    
    %     [atom(i).type]=atom(i).type{1}(1:2); % Elements do not have more than two characters;
    
    if strncmp(atom(i).type,{'OW'},2);atom(i).type={'Ow'};
    elseif strncmpi(atom(i).type,{'Hw'},2);atom(i).type={'Hw'};
    end
    
    if strncmpi(atom(i).type,{'Si'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'Sr'},2);atom(i).element={'Sr'};
    elseif strncmpi(atom(i).type,{'SY'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'SC'},2);atom(i).element={'Si'};
    elseif strncmp(atom(i).type,{'S'},1);atom(i).element={'S'};
    elseif strncmp(atom(i).type,{'st'},2);atom(i).element={'Si'};
    elseif strncmp(atom(i).type,{'s'},1);atom(i).element={'Si'};
    elseif strncmp(atom(i).type,{'Al'},2);atom(i).element={'Al'};
    elseif strncmp(atom(i).type,{'a'},1);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'AC'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'AY'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'Br'},2);atom(i).element={'Br'};
    elseif strncmpi(atom(i).type,{'B'},1);atom(i).element={'B'};
    elseif strncmpi(atom(i).type,{'I'},1);atom(i).element={'I'};
    elseif strncmpi(atom(i).type,{'Mg'},2);atom(i).element={'Mg'};
    elseif strncmpi(atom(i).type,{'Fe'},2);atom(i).element={'Fe'};
    elseif strncmpi(atom(i).type,{'F'},1);atom(i).element={'F'};
    elseif strncmpi(atom(i).type,{'U'},1);atom(i).element={'U'};
    elseif strncmpi(atom(i).type,{'V'},1);atom(i).element={'V'};
    elseif strncmpi(atom(i).type,{'Y'},1);atom(i).element={'Y'};
        %     elseif strncmpi(atom(i).type,{'Ow'},2);atom(i).element={water_O};
        %     elseif strncmpi(atom(i).type,{'Hw'},2);atom(i).element={water_H};
    elseif strncmpi(atom(i).type,{'Ow'},2);atom(i).element={'O'};
    elseif strncmpi(atom(i).type,{'Hw'},2);atom(i).element={'H'};
    elseif strncmpi(atom(i).type,{'Li'},2);atom(i).element={'Li'};
    elseif strncmpi(atom(i).type,{'Mn'},2);atom(i).element={'Mn'};
    elseif strncmpi(atom(i).type,{'Na'},2);atom(i).element={'Na'};
    elseif strncmpi(atom(i).type,{'Ni'},2);atom(i).element={'Ni'};
    elseif strncmpi(atom(i).type,{'Nh'},2);atom(i).element={'Nh'};
    elseif strncmpi(atom(i).type,{'Nb'},2);atom(i).element={'Nb'};
    elseif strncmpi(atom(i).type,{'Ne'},2);atom(i).element={'Ne'};
    elseif strncmpi(atom(i).type,{'No'},2);atom(i).element={'No'};
    elseif strncmpi(atom(i).type,{'N'},1);atom(i).element={'N'};
    elseif strncmpi(atom(i).type,{'K'},1);atom(i).element={'K'};
    elseif strncmpi(atom(i).type,{'Ca'},2);atom(i).element={'Ca'};
    elseif strncmp(atom(i).type,{'O'},1);atom(i).element={'O'};
    elseif strncmp(atom(i).type,{'o'},1);atom(i).element={'O'};
    elseif strncmp(atom(i).type,{'H'},1);atom(i).element={'H'};
    elseif strncmp(atom(i).type,{'h'},1);atom(i).element={'H'};
    elseif strncmpi(atom(i).type,{'Ca'},2);atom(i).element={'Ca'};
    elseif strncmpi(atom(i).type,{'Cl'},2);atom(i).element={'Cl'};
    elseif strncmpi(atom(i).type,{'C'},1);atom(i).element={'C'};
    else
        [atom(i).element{1}(1)]=upper(atom(i).type{1}(1));
        if size([atom(i).type{:}],2)>1
            [atom(i).element{1}(2)]=lower(atom(1).type{1}(2));
        end
    end
end


nAtoms=length(atom);
Atom_section=cell(nAtoms,10);

atom=orto_atom(atom,Box_dim);

%atom = number_type(atom);

% Write the file
fid = fopen(filename_out, 'wt');
fprintf(fid, '%s\r\n','data_matlab_gen');
fprintf(fid, '\r\n');
timestamp = datetime;
fprintf(fid, '%s     %s\r\n','_audit_creation_date',timestamp);
fprintf(fid, '%s\r\n','_audit_creation_method   generated by the Matlab ATOM scripts');

fprintf(fid, '\r\n');

fprintf(fid, '%-22s   %8.4f\r\n','_cell_length_a',a);
fprintf(fid, '%-22s   %8.4f\r\n','_cell_length_b',b);
fprintf(fid, '%-22s   %8.4f\r\n','_cell_length_c',c);
fprintf(fid, '%-22s   %8.4f\r\n','_cell_angle_alpha',alfa);
fprintf(fid, '%-22s   %8.4f\r\n','_cell_angle_beta',beta);
fprintf(fid, '%-22s   %8.4f\r\n','_cell_angle_gamma',gamma);
fprintf(fid, '%-22s   %8.4f\r\n','_cell_volume',Box_volume);

fprintf(fid, '\r\n');

fprintf(fid, '%s\r\n','_symmetry_cell_setting          triclinic');
fprintf(fid, '%s\r\n','_symmetry_space_group_name_Hall ''P 1''');
fprintf(fid, '%s\r\n','_symmetry_space_group_name_H-M  ''P 1''');
fprintf(fid, '%s\r\n','_symmetry_Int_Tables_number     1');

fprintf(fid, '\r\n');

fprintf(fid, '%s\r\n','_symmetry_equiv_pos_as_xyz ''x,y,z''');

fprintf(fid, '\r\n');

fprintf(fid, '%s\r\n','loop_');
%fprintf(fid, '%s\r\n','1 x,y,z');
% fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','loop_');
fprintf(fid, '%s\r\n','_atom_site_label');
fprintf(fid, '%s\r\n','_atom_site_type_symbol');
fprintf(fid, '%s\r\n','_atom_site_fract_x');
fprintf(fid, '%s\r\n','_atom_site_fract_y');
fprintf(fid, '%s\r\n','_atom_site_fract_z');

if isfield(atom,'charge')
    fprintf(fid, '%s\r\n','_atom_site_charge');
    for i = 1:nAtoms
        Atom_section(1:6) = [atom(i).type, atom(i).element, atom(i).xfrac, atom(i).yfrac, atom(i).zfrac, atom(i).charge];
        fprintf(fid,'%-11s%-6s%14.7f%14.7f%14.7f%12.5f\r\n',Atom_section{1:6});
    end
else
    for i = 1:nAtoms
        Atom_section(1:5) = [atom(i).type, atom(i).element, atom(i).xfrac, atom(i).yfrac, atom(i).zfrac];
        fprintf(fid,'%-11s%-6s%14.7f%14.7f%14.7f\r\n',Atom_section{1:5});
    end
end
fprintf(fid, '\r\n');
fprintf(fid, '\r\n');
% fprintf(fid, '%s\r\n','#END');

fclose(fid);

% for i = 1:nAtoms
%     Atom_section(1:13) = ['ATOM  ', atom(i).index, atom(i).type, atom(i).resname, 'A',atom(i).molid, atom(i).x, atom(i).y, atom(i).z,1,1,atom(i).element,atom(i).formalcharge];
%     fprintf(fid,'%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2i\n',Atom_section{1:13});
% end
% assignin('caller','Bond_index',Bond_index);
% assignin('caller','Angle_index',Angle_index);


