%% This script creates a gromacs .itp file from an .xyz file. If te .xyz
%% file contains two clay layers, this script writes two .itp files. It assumes
%% you have previously made a .xyz file from the clayff_system.m script
tic
clear all;
format compact;

file_title = 'Gromacs awesome itp file'; % Header in output file
molecule_name = 'MMT_4';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename='interface_MMT_4.gro';
filename_out=strcat(molecule_name,'_itp.gro');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot = 1; % 1/0 for yes/no
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import_atom(filename);
Atom_label = unique([atom.type]);
UCinX = 6; % Number of unit cells in the x direction
UCinY = 4; % Number of unit cells in the y direction
max_short_dist=1.25;      % Sets min bond distance
max_long_dist=2.3;       % Sets max bond distance
nrexcl=2; % See the gromacs manual
nAtoms = sum(strcmp([atom.resname],[atom(1).resname]));
precision = 7; % num2str(X,precision)
Natom_types = size(Atom_label,2);
Nbond_types = 2; % For MMT-O-H and H-O-H
Nangle_types = 1; % For H-O-H
INTERFACE_param(Atom_label,'SPC/E'); % Import forcefield parameters, state water model, SPC or SPC/E

%%

ind_H=find(strncmp([atom.type],{'H'},1));
ind_Oh=find(strncmp([atom.type],{'Oh'},2));
ind_Oct=[find(strncmp([atom.type],{'Al'},2)) find(strncmp([atom.type],{'Mgo'},3))];

fid = fopen(strcat(molecule_name,'.itp'), 'wt');

fprintf(fid, '%s % s\r\n',';',file_title);
fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','[ moleculetype ]');
fprintf(fid, '%s % s\r\n',';','molname   nrexcl');
fprintf(fid, '%s       %d\r\n',molecule_name,nrexcl);
fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','[ atoms ]');
fprintf(fid, '%s\r\n','; id   attype  resnr resname  atname   cgnr charge      mass');

atom=bond_angle_atom(atom,Box_dim,max_short_dist,max_long_dist);

Atom_label_ID=ones(size(atom,2),1);
for i = 1:size(XYZ_labels,1);
    if sum(ismember(Atom_label,[atom(i).type])) > 0;
        Atom_label_ID(i,1)=find(ismember(Atom_label,[atom(i).type])==1);
    end
    Atoms_data(i,:) = {i, char([atom(i).type]),1,molecule_name(1:3),char([atom(i).type]),i, Charge(Atom_label_ID(i,1)), Masses(Atom_label_ID(i,1))}; 
    fprintf(fid, '%-4i%6s%8i%8s%8s%8i\t%8.6f\t%8.6f\r\n', Atoms_data{i,:});
end

fprintf(fid, '\r\n');
fprintf(fid, '[ bonds ] \r\n');
fprintf(fid, '%s\r\n','; i    j    type');

count_b = 1;
bondtype=1; 
explicit_bonds = 1;
while count_b <= nBonds;
    if explicit_bonds == 1;
        if sum(ismember(Bond_index(count_b,1:2),ind_H))>0
            r=0.09290;
            kb=414216;
        else
            r=Bond_index(count_b,3)/10*1.05;
            kb=359824;
        end
        
        Bond_order(count_b,:)= {Bond_index(count_b,1), Bond_index(count_b,2), bondtype, r, kb, ';',strtrim(char(XYZ_labels(Bond_index(count_b,1)))), strtrim(char(XYZ_labels(Bond_index(count_b,2))))};
        fprintf(fid, '%-5i %-5i %-5i %-8.4f %-8.4f %s %s-%s\r\n', Bond_order{count_b,:});
        count_b = count_b + 1;
    else
        Bond_order(count_b,:)= {Bond_index(count_b,1), Bond_index(count_b,2), bondtype, ';',strtrim(char(XYZ_labels(Bond_index(count_b,1)))), strtrim(char(XYZ_labels(Bond_index(count_b,2)))), Bond_index(count_b,3)/10};
        fprintf(fid, '%-5i %-5i %-5i %s %s-%s %-8.4f\r\n', Bond_order{count_b,:});
        count_b = count_b + 1;
    end
end

disp('These atom types has bonds')
unique(Bond_order(:,7:8))

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '[ angles ] \r\n');
fprintf(fid, '%s\r\n','; i    j   k   type');

%     %% To remove angles with 'Al'
%     Al_ind=find(strcmp(XYZ_labels(:,1),'Al'));
%     [Al_row,Al_col]=ind2sub(size(Angle_index),find(ismember(Angle_index,Al_ind)));
%     Angle_index(Al_row,:)=[];
%     %% To remove angles with 'Si'
%     Si_ind=find(strcmp(XYZ_labels(:,1),'Si'));
%     [Si_row,Si_col]=ind2sub(size(Angle_index),find(ismember(Angle_index(:,2),Si_ind)));
%     Angle_index(Si_row,:)=[];

count_a = 1;explicit_angles = 1;
angletype=1; Angle_order={};
%     Angle_index=sortrows(Angle_index); %% Why???
while count_a <= length(Angle_index); %nAngles;
    if explicit_angles == 1;
        
        if sum(ismember(Angle_index(count_a,1:3),ind_H))>0 && sum(ismember(Angle_index(count_a,1:3),ind_Oh))>0 && sum(ismember(Angle_index(count_a,1:3),ind_Oct))>0
            adeg=126; % From most recent CHARMM prm file 116.2;
            ktheta=376.56; % since 45*4.184*2;% earlier 96.232*10; %
        else
            adeg=Angle_index(count_a,4);
            ktheta=1422.56;
        end
        
        Angle_order(count_a,:)= {Angle_index(count_a,1), Angle_index(count_a,2), Angle_index(count_a,3), angletype, adeg,	ktheta, ';', strtrim(char(XYZ_labels(Angle_index(count_a,1)))), strtrim(char(XYZ_labels(Angle_index(count_a,2)))), strtrim(char(XYZ_labels(Angle_index(count_a,3))))};
        fprintf(fid, '%-5i %-5i %-5i %-5i %-8.4f %-8.4f %s %s-%s-%s\r\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    else
        Angle_order(count_a,:)= {Angle_index(count_a,1), Angle_index(count_a,2), Angle_index(count_a,3), angletype, ';', strtrim(char(XYZ_labels(Angle_index(count_a,1)))), strtrim(char(XYZ_labels(Angle_index(count_a,2)))), strtrim(char(XYZ_labels(Angle_index(count_a,3)))), Angle_index(count_a,4)};
        fprintf(fid, '%-5i %-5i %-5i %-5i %s %s-%s-%s %-8.4f\r\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    end
end

%% Defining [ exclusions ]
%     if length(Angle_index) > 0;
%
%         fprintf(fid, '\r\n');
%         fprintf(fid, '\r\n');
%         fprintf(fid, '[ exclusions ] \r\n');
%         fprintf(fid, '%s\r\n','; i    j   k   type');
%
%         count_excl = 1;
%         Excl_index=[Angle_index(:,1) Angle_index(:,2) Angle_index(:,3); Angle_index(:,2) Angle_index(:,3) Angle_index(:,1); Angle_index(:,2) Angle_index(:,3) Angle_index(:,1)];
%         while count_excl <= length(Excl_index);
%             Excl_order(count_excl,:)= {Excl_index(count_excl,1), Excl_index(count_excl,2), Excl_index(count_excl,3),';', strtrim(char(XYZ_labels(Excl_index(count_excl,1)))), strtrim(char(XYZ_labels(Excl_index(count_excl,2)))), strtrim(char(XYZ_labels(Excl_index(count_excl,3))))};
%             fprintf(fid, '%-5i %-5i %-5i %s %s-%s-%s\r\n', Excl_order{count_excl,:});
%             count_excl = count_excl + 1;
%         end
%
%     end

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '#ifdef POSRES_Y_10 \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:40:nAtoms;
    pos_res(i,:) = {num2str(i), '1', '1000', '10', '1000'};
    fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
    fprintf(fid, '\r\n');
end
fprintf(fid, '#endif \r\n');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '#ifdef POSRES_Y_100 \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:40:nAtoms;
    pos_res(i,:) = {num2str(i), '1', '1000', '100', '1000'};
    fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
    fprintf(fid, '\r\n');
end
fprintf(fid, '#endif \r\n');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '#ifdef POSRES_INTERFACE \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms;
    pos_res(i,:) = {num2str(i), '1', '1000', '1000', '1000'};
    fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
    fprintf(fid, '\r\n');
end
fprintf(fid, '#endif \r\n');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '#ifdef POSRES_noH \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms;
    if strncmp(XYZ_labels(i),'H',1)==0;
        pos_res(i,:) = {num2str(i), '1', '1000', '1000', '1000'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
        fprintf(fid, '\r\n');
    end
end
fprintf(fid, '#endif \r\n');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '#ifdef POSRES_XYZ \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:40:nAtoms;
    pos_res(i,:) = {num2str(i), '1', '1000', '1000', '1000'};
    fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
    fprintf(fid, '\r\n');
end
fprintf(fid, '#endif \r\n');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '#ifdef POSRES_XY \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:40:nAtoms;
    pos_res(i,:) = {num2str(i), '1', '1000', '1000', '0'};
    fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
    fprintf(fid, '\r\n');
end
fprintf(fid, '#endif \r\n');

fclose(fid);

Total_charge = check_INTERFACE_charge(atom)
%% check these out!!!
dist_matrix = dist_matrix_atom(atom,Box_dim);
wrapped_atom = analyze_atom(atom,Box_dim,1.25,2.15);

% [atom(strcmp([atom.type],{'Ow'})).type]=deal({'OW'});
% [atom(strcmp([atom.type],{'Hw'})).type]=deal({'HW'});
atom=center_atom(atom,Box_dim,'all','xyz');
atom=translate_atom(atom,[0 0 -Box_dim(3)/2],'all')
%%%%%% Write structure to file %%%%%%
% write_atom_gro(atom,Box_dim,filename_out);

%%%%%% Plot the structure %%%%%%%%%%%
if plot==1;vmd(atom,Box_dim);end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


toc
