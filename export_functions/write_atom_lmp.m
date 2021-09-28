%% write_atom_lmp.m
% * This script creates and prints a lammps data file (.lj). Works best for
% Clayff systems, defining the ffname 'clayff' and watermodel 'spce'.
% Nevertheless, this new version should be able to handle bonds|angles|dihedrals
%
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_lmp(atom,Box_dim,filename,1.25,1.25,'clayff','spce') % Basic input arguments
%
function write_atom_lmp(atom,Box_dim,filename,varargin)
format compact;

% Change these values in case you want to merge a lammps topology file with
% another one already haveing bonds, angles, dihedrals...
prev_atom_index=0;
prev_atom_types=0;
prev_bond_num=0;
prev_bond_types=0;
prev_angle_num=0;
prev_angle_types=0;
prev_dihedral_num=0;
prev_dihedral_types=0;
prev_mol_index=0;

precision = 7; % num2str(X,precision)
ind_HW=strncmpi([atom.type],'Ow',2);[atom(ind_HW).type]=deal({'Ow'});
ind_HW=strncmpi([atom.type],'Hw',2);[atom(ind_HW).type]=deal({'Hw'});
Atom_labels=sort(unique([atom.type]));
natom_labels=size(Atom_labels,2);

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
    if strcmpi(ffname,'clayff_2004')
        clayff_2004_param(sort(unique([atom.type])),watermodel);
        if ~isfield(atom,'charge')
            atom = charge_atom(atom,Box_dim,'clayff_2004',watermodel,'adjust');
        end
        Total_charge=sum([atom.charge])
        round(Total_charge,5)
    elseif strncmpi(ffname,'clayff',5)
        clayff_param(sort(unique([atom.type])),watermodel);
        atom = charge_atom(atom,Box_dim,'clayff',watermodel,'more');
        Total_charge
        round(Total_charge,5)
    end
else
    %     disp('Unknown forcefield, will try clayff')
    %     clayff_param(sort(unique([atom.type])),'spc/e');
    %     atom = charge_atom(atom,Box_dim,'clayff','spc/e');
    %     Total_charge
    disp('Forcefield not stated, will make some assumptions then...')
    pause(2)
    ffname='clayff_2004'
    watermodel='SPC/E'
    clayff_2004_param(sort(unique([atom.type])),watermodel);
    if ~isfield(atom,'charge')
        atom = charge_atom(atom,Box_dim,'clayff_2004',watermodel,'adjust');
    end
    Total_charge=sum([atom.charge])
    round(Total_charge,5)
end

% Scan the xyz data and look for O-H bonds and angles
atom=bond_angle_dihedral_atom(atom,Box_dim,short_r,long_r);
assignin('caller','atom',atom);
assignin('caller','Bond_index',Bond_index);
assignin('caller','Angle_index',Angle_index);
assignin('caller','Dihedral_index',Dihedral_index);
nAtoms=size(atom,2)
nBonds
nAngles
nDihedrals
 
if nBonds>0
    % Calculate the bond_types
    bond_pairs=[[atom(Bond_index(:,1)).type]' [atom(Bond_index(:,2)).type]'];
    % b1=join([bond_pairs(:,1) bond_pairs(:,2)]); % Does not work in older MATLAB versions?
    b1=strcat(bond_pairs(:,1),{' '},bond_pairs(:,2));
    b1=unique(b1,'stable')
    % b2=join([bond_pairs(:,2) bond_pairs(:,1)]); % Does not work in older MATLAB versions?
    b2=strcat(bond_pairs(:,2),{' '},bond_pairs(:,1));
    b2=unique(b2,'stable');
    nbond_types=size(b1,1);
    bond_pairs=join(bond_pairs);
    for i=1:size(bond_pairs,1)
        [ind,bond_types(i)]=ismember(bond_pairs(i),[b1;b2]);
    end
    bond_types(bond_types>nbond_types)=bond_types(bond_types>nbond_types)-nbond_types;
else
    nbond_types=0;
end

if nAngles>0
    % Calculate the angle_types
    angle_triplets=[[atom(Angle_index(:,1)).type]' [atom(Angle_index(:,2)).type]' [atom(Angle_index(:,3)).type]'];
    a1=join([angle_triplets(:,1) angle_triplets(:,2) angle_triplets(:,3)]);a1=unique(a1,'stable');
    a2=join([angle_triplets(:,3) angle_triplets(:,2) angle_triplets(:,1)]);a2=unique(a2,'stable');
    nangle_types=size(a1,1);
    angle_triplets=join(angle_triplets);
    for i=1:size(angle_triplets,1)
        [ind,angle_types(i)]=ismember(angle_triplets(i),[a1;a2]);
    end
    angle_types(angle_types>nangle_types)=angle_types(angle_types>nangle_types)-nangle_types;
else
    nangle_types=0;
end

if nDihedrals>0
    % Calculate the dihedral_types
    dihedral_quads=[[atom(Dihedral_index(:,1)).type]' [atom(Dihedral_index(:,2)).type]' [atom(Dihedral_index(:,3)).type]' [atom(Dihedral_index(:,4)).type]'];
    d1=join([dihedral_quads(:,1) dihedral_quads(:,2) dihedral_quads(:,3) dihedral_quads(:,4)]);d1=unique(d1,'stable');
    d2=join([dihedral_quads(:,4) dihedral_quads(:,3) dihedral_quads(:,2) dihedral_quads(:,1)]);d2=unique(d2,'stable');
    ndihedral_types=size(d1,1);
    dihedral_quads=join(dihedral_quads);
    for i=1:size(dihedral_quads,1)
        [ind,dihedral_types(i)]=ismember(dihedral_quads(i),[d1;d2]);
    end
    dihedral_types(dihedral_types>ndihedral_types)=dihedral_types(dihedral_types>ndihedral_types)-ndihedral_types;
else
    ndihedral_types=0;
end

% End of settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('Box_dim','var');disp('We need to set Box_dim?');
    % Box_dim = ?
    pause
end

lx=Box_dim(1);ly=Box_dim(2);lz=Box_dim(3);
if length(Box_dim)>3
    triclinic = 1; % 1 or 0 for ortoghonal
    xy=Box_dim(6);xz=Box_dim(8);yz=Box_dim(9);
else
    triclinic = 0;
    xy=0;xz=0;yz=0;
end


% Some settings needed in order to print a proper lammps in file %%

% Start printing the lammps .lj file
fid = fopen(filename, 'wt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_title = strcat('LAMMPS input data file #',datestr(now)); % Header in output file
fprintf(fid, '%s\n', file_title);
fprintf(fid, '\n');

%%
AtomBondAnglestring = {num2str(size(atom,2)), 'atoms';...
    num2str(nBonds), 'bonds';...
    num2str(nAngles), 'angles';...
    num2str(nDihedrals), 'dihedrals';...
    ' ',' ';...
    num2str(natom_labels), 'atom types';...
    num2str(nbond_types), 'bond types';...
    num2str(nangle_types), 'angle types';...
    num2str(ndihedral_types), 'dihedral types'};

Boxsizestring = {  num2str(0,precision), num2str(lx,precision), 'xlo', 'xhi';...
    num2str(0,precision), num2str(ly,precision), 'ylo', 'yhi';...
    num2str(0,precision), num2str(lz,precision), 'zlo', 'zhi';...
    num2str(xy,precision), num2str(xz,precision),num2str(yz,precision), 'xy xz yz'};

for i = 1:size(AtomBondAnglestring,1)
    fprintf(fid, '%-s %-s\n', AtomBondAnglestring{i,:});
end
fprintf(fid, '\n');

for i = 1:size(Boxsizestring,1) % Im not sure what this code is doing
    if triclinic == 0
        Boxsizestring(end,:) = {' ',' ',' ',' '};
    end
    fprintf(fid, '%-s %-s %-s %-s\n', Boxsizestring{i,:});
end
%%
% fprintf(fid, '\n');
fprintf(fid, 'Masses \n');
fprintf(fid, '\n');

% if isfield(atom,'mass')
%     for i =1:natom_labels
%         masses(i,:) = {i+prev_atom_types, num2str([atom(i).mass],precision), '#', Atom_labels{i}};
%         fprintf(fid, '%i %-s %s %s\n', masses{i,:});
%     end
% else exist('Masses','var');
    for i =1:natom_labels
        masses(i,:) = {i+prev_atom_types, num2str(Masses(i),precision), '#', Atom_labels{i}};
        fprintf(fid, '%i %-s %s %s\n', masses{i,:});
    end
% end

fprintf(fid, '\n');
%%
%%
fprintf(fid, 'Pair Coeffs \n');
fprintf(fid, '\n');

for i =1:natom_labels
    paircoeffs(i,:) = {i, Epsilon(i)*kcalmol_to_eV, Sigma(i), '#', Atom_labels{i}};
    fprintf(fid, '%i %E %f %s %s\n', paircoeffs{i,:});
end
fprintf(fid, '\n');
%
fprintf(fid, 'Atoms \n');
fprintf(fid, '\n');

for i = 1:nAtoms
    if sum(ismember(Atom_labels,[atom(i).type])) > 0
        Atom_label_ID(i,1)=find(ismember(Atom_labels,[atom(i).type])==1);
    end
    Atoms_data(i,:) = {i+prev_atom_index, [atom(i).molid]+prev_mol_index, Atom_label_ID(i,1), [atom(i).charge],[atom(i).x],[atom(i).y],[atom(i).z]};
    fprintf(fid, '\t%-i\t%-i\t%-i\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n', Atoms_data{i,:});
end

%%
fprintf(fid, '\n');
fprintf(fid, '\n');

if nBonds>0
    % Prints bond data
    fprintf(fid, 'Bonds \n');
    fprintf(fid, '\n');
    count_b = 1;
    while count_b <= nBonds
        Bond_order(count_b,:)= {count_b+prev_bond_num, bond_types(count_b)+prev_bond_types, Bond_index(count_b,1)+prev_atom_index, Bond_index(count_b,2)+prev_atom_index};
        fprintf(fid, '\t%-i %-i %-i %-i\n', Bond_order{count_b,:});
        count_b = count_b + 1;
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');
end

if nAngles>0
    % Prints angle data
    fprintf(fid, 'Angles \n');
    fprintf(fid, '\n');
    count_a = 1;
    while count_a <= nAngles
        Angle_order(count_a,:)= {count_a+prev_angle_num, angle_types(count_a)+prev_angle_types, Angle_index(count_a,1)+prev_atom_index,Angle_index(count_a,2)+prev_atom_index,Angle_index(count_a,3)+prev_atom_index};
        fprintf(fid, '\t%-i %-i %-i %-i %-i\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');
end

if nDihedrals>0
    % Prints dihedral data
    fprintf(fid, 'Dihedrals \n');
    fprintf(fid, '\n');
    count_d = 1;
    while count_d <= nDihedrals
        Dihedral_order(count_d,:)= {count_d+prev_dihedral_num, dihedral_types(count_d)+prev_dihedral_types, Dihedral_index(count_d,1)+prev_atom_index,Dihedral_index(count_d,2)+prev_atom_index,Dihedral_index(count_d,3)+prev_atom_index,Dihedral_index(count_d,4)+prev_atom_index};
        fprintf(fid, '\t%-i %-i %-i %-i %-i %-i\n', Dihedral_order{count_d,:});
        count_d = count_d + 1;
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');
end

fclose(fid);

