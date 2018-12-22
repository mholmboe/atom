%% write_atom_lmp.m
% * This script creates and prints a lammps data file (.lj). Works best for
% Clayff systems
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_lmp(atom,Box_dim,filename,1.25,1.25,'clayff','spce')
%
function write_atom_lmp(atom,Box_dim,filename,varargin)
format long;
prev_atom_index=0;
prev_atom_types=0;
prev_bond_num=0;
prev_bond_types=0;
prev_angle_num=0;
prev_angle_types=0;
prev_mol_index=0;

nAtoms=size(atom,2);

if regexp(filename,'.lj') ~= false
    filename = filename;
else
    filename = strcat(filename,'.lj');
end

if nargin > 3
    short_r=cell2mat(varargin(1));
    long_r=cell2mat(varargin(2));
else
    short_r=1.25;
    long_r=1.25; % long=short since clayff
end

if nargin>5
    ffname=varargin(3);
    if nargin>6
        watermodel=varargin(4)
    else
        disp('Unknown watermodel, will try SPC/E')
        watermodel='SPC/E';
    end
    if strncmpi(ffname,'clayff',5)
        clayff_param(sort(unique([atom.type])),watermodel);
        atom = charge_atom(atom,Box_dim,'clayff',watermodel);
        Total_charge
    end
else
    disp('Unknown forcefield, will try clayff')
    clayff_param(sort(unique([atom.type])),'spc/e');
    atom = charge_atom(atom,Box_dim,'clayff','spc/e');
    Total_charge
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Title=filename;
file_title = strcat('LAMMPS input data file #',datestr(now)); % Header in output file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some settings needed in order to print a proper lammps in file %%
bond_atoms={};angle_atoms={};
if sum(strcmpi([atom.type],'Oh')); bond_atoms={'Oh'}; end
if sum(strcmpi([atom.type],'Ohmg')); bond_atoms={bond_atoms{:} 'Ohmg'}; end
if sum(strcmpi([atom.type],'OW')); bond_atoms={bond_atoms{:} 'OW'}; end
if sum(strcmpi([atom.type],'OW')); angle_atoms={'OW'};end
% End of settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('Box_dim','var');disp('We need to set Box_dim?');
    % Box_dim = ?
    pause
end

Natoms = size(atom,2);
precision = 7; % num2str(X,precision)
ind_HW=strncmpi([atom.type],'Ow',2);[atom(ind_HW).type]=deal({'Ow'})
ind_HW=strncmpi([atom.type],'Hw',2);[atom(ind_HW).type]=deal({'Hw'})
Atom_label=sort(unique([atom.type]));
Natom_types = length(Atom_label);
if strncmpi(ffname,'clayff',5)
    clayff_param(Atom_label,watermodel); % Import forcefield parameters, state water model, SPC or SPC/E
else
    disp('Only clayff implemented sofar')
end
lx=Box_dim(1);ly=Box_dim(2);lz=Box_dim(3);
if length(Box_dim)>3
    triclinic = 1; % 1 or 0 for ortoghonal
    xy=Box_dim(6);xz=Box_dim(8);yz=Box_dim(9);
else
    triclinic = 0;
    xy=0;xz=0;yz=0;
end

% Scan the xyz data and look for O-H bonds and angles
atom=bond_angle_atom(atom,Box_dim,short_r,long_r);

% Start printing the lammps .lj file
fid = fopen(filename, 'wt');

fprintf(fid, '%s\r\n', file_title);
fprintf(fid, '\r\n');

%%
AtomBondAnglestring = {num2str(size(atom,2)), 'atoms';...
    num2str(size(Bond_index,1)), 'bonds';... %num2str(Nbonds), 'bonds';...
    num2str(size(Angle_index,1)), 'angles';...%num2str(Nangles), 'angles';...
    ' ',' ';...
    num2str(length(Atom_label)), 'atom types';...
    num2str(length(bond_atoms)), 'bond types';...
    num2str(length(angle_atoms)), 'angle types'};

Boxsizestring = {  num2str(0,precision), num2str(lx,precision), 'xlo', 'xhi';...
    num2str(0,precision), num2str(ly,precision), 'ylo', 'yhi';...
    num2str(0,precision), num2str(lz,precision), 'zlo', 'zhi';...
    num2str(xy,precision), num2str(xz,precision),num2str(yz,precision), 'xy xz yz'};

for i = 1:size(AtomBondAnglestring,1)
    fprintf(fid, '%-s %-s\r\n', AtomBondAnglestring{i,:});
end
fprintf(fid, '\r\n');

for i = 1:size(Boxsizestring,1) % Im not sure what this code is doing
    if triclinic == 0
        Boxsizestring(end,:) = {' ',' ',' ',' '};
    end
    fprintf(fid, '%-s %-s %-s %-s\r\n', Boxsizestring{i,:});
end
%%
fprintf(fid, '\r\n');
fprintf(fid, 'Masses \r\n');
fprintf(fid, '\r\n');

for i =1:length(Atom_label)
    masses(i,:) = {i+prev_atom_types, num2str(Masses(i),precision), '#', Atom_label{i}};
    fprintf(fid, '%i %-s %s %s\r\n', masses{i,:});
end
fprintf(fid, '\r\n');
%%
%%
fprintf(fid, 'Pair Coeffs \r\n');
fprintf(fid, '\r\n');

for i =1:length(Atom_label)
    paircoeffs(i,:) = {i, Epsilon(i)*kcalmol_to_eV, Sigma(i), '#', Atom_label{i}};
    fprintf(fid, '%i %E %f %s %s\r\n', paircoeffs{i,:});
end
fprintf(fid, '\r\n');
%
fprintf(fid, 'Atoms \r\n');
fprintf(fid, '\r\n');
% Construct and print the atoms properties section
% lmp_atom_style_full_func(fid,Atom_labels,Charge,[atom.type]',[[atom.x]' [atom.y]' [atom.z]']);

% atom = bond_angle_atom(atom,Box_dim,max_short_dist,max_long_dist);
% atom = charge_clayff_atom(atom,Box_dim); %{'Al' 'Mgo' 'Si' 'H'},(1.575 1.36 2.1 0.425));

atom = charge_atom(atom,Box_dim,'clayff',watermodel);

for i = 1:size(atom,2)
    %     Atoms_data(i,:) = {i, molID(i), Atom_label_ID(i), Charge(Atom_label_ID(i,1)), XYZ_data(i,1),XYZ_data(i,2), XYZ_data(i,3)};
    if sum(ismember(Atom_label,[atom(i).type])) > 0
        Atom_label_ID(i,1)=find(ismember(Atom_label,[atom(i).type])==1);
    end
    Atoms_data(i,:) = {i+prev_atom_index, [atom(i).molid]+prev_mol_index, Atom_label_ID(i,1), [atom(i).charge],[atom(i).x],[atom(i).y],[atom(i).z]};
    fprintf(fid, '\t%-i\t%-i\t%-i\t%8.5f\t%8.5f\t%8.5f\t%8.5f\r\n', Atoms_data{i,:});
end

% Atom_prop = {atomID(1:end-1), molID(1:end-1), Atom_label_ID(:,1), Charge(1,:), XYZ_data(i,1),XYZ_data(i,2), XYZ_data(i,3)};
% Atom_label_ID=ones(size(atom,2),1);
% for i = 1:nAtoms;
%     if sum(ismember(Atom_label,[atom(i).type])) > 0;
%         Atom_label_ID(i,1)=find(ismember(Atom_label,[atom(i).type])==1);
%     end
%     Atoms_data(i,:) = {i, char([atom(i).type]),1,molecule_name(1:3),char([atom(i).type]),i, [atom(i).charge], Masses(Atom_label_ID(i,1))};
%     fprintf(fid, '%-4i%6s%8i%8s%8s%8i\t%8.5f\t%8.6f\r\n', Atoms_data{i,:});
% end

%%
fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

% Prints bond data
fprintf(fid, 'Bonds \r\n');
fprintf(fid, '\r\n');

count_b = 1;
while count_b <= nBonds
    %     if Bond_index(count_b,1) <= solute_atoms*layers;
    %         bondtype=1;
    %     else
    %         bondtype=2;
    %     end
    for i=1:length(bond_atoms)
        if strcmpi([atom(Bond_index(count_b,1)).type],bond_atoms(i))
            bondtype=i+prev_bond_types;
        else
            bondtype=1+prev_bond_types;
        end
    end
    Bond_order(count_b,:)= {count_b+prev_bond_num, bondtype+prev_bond_types, Bond_index(count_b,1)+prev_atom_index, Bond_index(count_b,2)+prev_atom_index};
    fprintf(fid, '\t%-i %-i %-i %-i\r\n', Bond_order{count_b,:});
    count_b = count_b + 1;
end

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

% Prints angle data
fprintf(fid, 'Angles \r\n');
fprintf(fid, '\r\n');

if length(angle_atoms)>0
    count_a = 1;
    while count_a <= nAngles
        %     angletype = 1;
        for i=1:length(angle_atoms)
            if strcmpi([atom(Angle_index(count_a,1)).type],angle_atoms(i))
                angletype=i;
            else
                angletype=1;
            end
        end
        Angle_order(count_a,:)= {count_a+prev_angle_num, angletype+prev_angle_types, Angle_index(count_a,1)+prev_atom_index,Angle_index(count_a,2)+prev_atom_index,Angle_index(count_a,3)+prev_atom_index};
        fprintf(fid, '\t%-i %-i %-i %-i %-i\r\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    end
    fprintf(fid, '\r\n');
    fprintf(fid, '\r\n');
end

fclose(fid);

