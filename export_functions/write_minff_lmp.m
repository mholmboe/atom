%% write_minff_lmp.m
% * This script creates and prints a lammps data file (.data). Works best for
% MINFF systems. Nevertheless, this new version should be able to handle bonds|angles|dihedrals
%
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_lmp(atom,Box_dim,filename) % Basic input arguments
% # write_atom_lmp(atom,Box_dim,filename,ff) % ff being a MATLAB struct containing forcefield parameters
%

function write_minff_lmp(atom,Box_dim,filename,varargin)
format compact;

% Some settings needed in order to print a proper lammps in file %%

bHdist=0.9572;              % Å
kbH=441050 /(2*4.184*10^2); % Gromacs Energy [kJ/mol] and distance [nm] to Lammps real units [kcal/mol] and [Å]  / 2 / 4.184 / 100.
kbM=0 /(2*4.184*10^2);      % Gromacs Energy [kJ/mol] and distance [nm] to Lammps real units [kcal/mol] and [Å]  / 2 / 4.184 / 100.
KANGLEH=125.52 /2/4.184;    % Gromacs [kJ/mol] to Lammps real units [kcal/mol]  / 2 / 4.184. Note the factor of 2. From Pouvrea, Greathouse, Cygan, and Andrey G. Kalinichev 2017
KANGLE=500.00 /2/4.184;     % Gromacs [kJ/mol] to Lammps real units [kcal/mol]  / 2 / 4.184. Note the factor of 2
angleH=110;                 % From Pouvrea, Greathouse, Cygan, and Andrey G. Kalinichev 2017
precision = 7;              % num2str(X,precision)

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

ind_OW=strncmpi([atom.type],'Ow',2);[atom(ind_OW).type]=deal({'Ow'});
ind_HW=strncmpi([atom.type],'Hw',2);[atom(ind_HW).type]=deal({'Hw'});
Atom_labels=sort(unique([atom.type]));
natom_labels=size(Atom_labels,2);

if nargin>3
    ff=varargin{1};
    if ~isstruct(ff)
        if regexp(ff,'.mat') ~= false
            ff = ff;
        else
            ff = strcat(ff,'.mat');
        end
        load(ff);
    end
    if size(ff,2)==1
        ff=ff.ff;
    end
    for i=1:size(Atom_labels,2)
        ind=strcmpi([ff.type],Atom_labels(i));
        ind
        Atom_labels(i)
        ff
        Epsilon(i)=ff(ind).e_kcalmol; % kcal/mol
        Sigma(i)=ff(ind).sigma_A; % Ångstrom
    end
else
    ff=[];
end

if nargin>4
    maxrshort=varargin{2};
    maxrlong=varargin{3};
else
    maxrshort=1.25;
    maxrlong=2.45;
end

if regexp(filename,'.data') ~= false
    filename = filename;
else
    filename = strcat(filename,'.data');
end

ffname='minff';
watermodel='OPC3'; % SPC/E, depreceateda
atom = mass_atom(atom);
Masses=[];
for i=1:size(Atom_labels,2)
    ind=find(strcmpi([atom.type],Atom_labels(i)));
    Masses(i)=[atom(ind(i)).mass];
end


if ~isfield(atom,'charge')
    atom=charge_minff_atom(atom,Box_dim,{'Al' 'Alt' 'Ale' 'Tio' 'Feo' 'Fet' 'Fee' 'Fe2' 'Fe2e' 'Fe3e' 'Na' 'K' 'Cs' 'Mgo' 'Mgh' 'Mge' 'Cao' 'Cah' 'Sit' 'Si' 'Sio' 'Site' 'Lio' 'H'},[1.782 1.782 1.985 2.48 1.14 1.14 1.14 0.7 0.86666 1.45 1 1 1 1.562 1.74 1.635 1.66 1.52 1.884 1.884 1.884 2.413 0.86 0.4]);
end

if exist('Total_charge','var')
    disp('Total charge for the .itp file was')
    round(Total_charge,5)
end

%% Find atomtype specific indexes
Bond_index=[];Angle_index=[];Dihedral_index=[];

ind_Hneighbours = find(~cellfun(@isempty,regexpi([atom.type],'h')));
ind_H=find(strncmpi([atom.type],{'H'},1));
ind_O=find(strncmpi([atom.type],{'O'},1));
ind_Osih=find(strncmpi([atom.type],{'Osih'},4));
ind_Alhh=find(strncmpi([atom.type],{'Oalhh'},5));
ind_Mghh=find(strncmpi([atom.type],{'Omhh'},4));
ind_Fehh=find(strncmpi([atom.type],{'Ofehh'},5));
ind_Oh=intersect(ind_O,ind_Hneighbours);
ind_Al=find(strncmpi([atom.type],'Al',2));
ind_Al=find(strcmp([atom.type],'Al'));
ind_Mgo=find(ismember([atom.type],{'Mgo' 'Mgh'}));
ind_Si=find(strncmpi([atom.type],{'Si'},2));
ind_Oct=sort([ind_Al ind_Mgo]);
ind_Edge=unique([ind_H ind_Alhh ind_Mghh ind_Fehh ind_Osih]);

atom = bond_angle_atom(atom,Box_dim,maxrshort,maxrlong);

%     %% To only keep bonds to atoms also bonded to H's, uncomment the next four lines
%     disp('Keeping only bonds with H')
%     [h_row,h_col]=ind2sub(size(Bond_index),find(ismember(Bond_index,ind_Hneighbours)));
%     Bond_index=Bond_index(h_row,:);
%     nBonds=size(Bond_index,1);

%% To only keep bonds to H's, uncomment the next three lines
%     [H_row,H_col]=ind2sub(size(Bond_index),find(ismember(Bond_index,ind_H)));
%     Bond_index=Bond_index(H_row,:);
%     nBonds=size(Bond_index,1);

%     %% To only keep edge bonds (and all O-H), uncomment the next four lines
disp('Keeping only bonds with H or edge-O ')
[h_row,h_col]=ind2sub(size(Bond_index),find(ismember(Bond_index,ind_Edge)));
Bond_index=Bond_index(h_row,:);

%% To only keep bonds between Osih - H, uncomment the next four lines
%     disp('Keeping only bonds with H')
%     [h_row,h_col]=ind2sub(size(Bond_index),find(ismember(Bond_index,ind_Osih)));
%     Bond_index=Bond_index(h_row,:);
%     nBonds=size(Bond_index,1);

%     %% To remove bonds with 'Al'
%     [Al_row,Al_col]=ind2sub(size(Bond_index),find(ismember(Bond_index,ind_Al)));
%     Bond_index(Al_row,:)=[];
%     nBonds=size(Bond_index,1);

%     %% To remove bonds with 'Si'
%     [Si_row,Si_col]=ind2sub(size(Bond_index),find(ismember(Bond_index(:,2),ind_Si)));
%     Bond_index(Si_row,:)=[];
%     nBonds=size(Bond_index,1);

%    %% To remove bonds larger than certain rmin, uncomment next two lines
%     rm_ind=find(Bond_index(:,3)>1.25);
%     Bond_index(rm_ind,:)=[];
%     nBonds=size(Bond_index,1);

[Y,I]=sort(Bond_index(:,1));
Bond_index=Bond_index(I,:);
Bond_index = unique(Bond_index,'rows','stable');

% if strncmpi(ffname,'minff',5)
disp('Keeping all angles with O... ')
%     disp('What to do with edge angles with Al')
%     disp('Removing angles with Al')
%     [Al_row,Al_col]=ind2sub(size(Angle_index),find(ismember(Angle_index,ind_Al)));
%     Angle_index(Al_row,:)=[];
%     % To remove angles with 'Si'
%     Si_ind=find(strcmp(XYZ_labels(:,1),'Si'));
%     [Si_row,Si_col]=ind2sub(size(Angle_index),find(ismember(Angle_index(:,2),Si_ind)));
%     Angle_index(Si_row,:)=[];
% %% Extra stuff
%     Angle_index
%     rm_ind=find(Angle_index(:,4)>150|Angle_index(:,4)<60);
%     Angle_index(rm_ind,:)=[];
%     [Y,I]=sort(Angle_index(:,2));
%     Angle_index=Angle_index(I,:);
%     Angle_index = unique(Angle_index,'rows','stable');
% end

assignin('caller','atom',atom);
% assignin('caller','Epsilon',Epsilon);
% assignin('caller','Sigma',Sigma);
assignin('caller','Masses',Masses);
assignin('caller','Bond_index',Bond_index);
assignin('caller','Angle_index',Angle_index);
assignin('caller','Dihedral_index',Dihedral_index);
nAtoms=size(atom,2)
nBonds=size(Bond_index,1)
nAngles=size(Angle_index,1)
nDihedrals=size(Dihedral_index,1)

if nBonds>0

    bond_pairs=[[atom(Bond_index(:,1)).type]' [atom(Bond_index(:,2)).type]'];
    bond_info1=[];%bond_info2=[];
    for i = 1:size(bond_pairs,1)
        bond_info1{i,1} = sprintf('%s %s %.3f', bond_pairs{i,1}, bond_pairs{i,2}, Bond_index(i,3));
        % bond_info2{i,1} = sprintf('%s %s %.3f', bond_pairs{i,2}, bond_pairs{i,1}, Bond_index(i,3));
    end
    b1=bond_info1;[b1,idxb1]=unique(b1,'stable');
    % b2=bond_info1;b2=unique(b2,'stable');
    nbond_types=size(b1,1);
    for i=1:size(bond_info1,1)
        [ind,bond_types(i)]=ismember(bond_info1(i),[b1]); %;b2]);
        % bond_coeffs
    end
    bond_types(bond_types>nbond_types)=bond_types(bond_types>nbond_types)-nbond_types;

    bond_coeffs=[1:length(idxb1)]';
    bond_dist=Bond_index(idxb1,3);
    bond_kb=Bond_index(idxb1,3);
    bond_kb(bond_dist>1.25)=kbM;
    bond_kb(bond_dist<1.25)=kbH;

    bond_dist(bond_dist<1.25)=bHdist;
    bond_coeffs=[bond_coeffs bond_kb bond_dist];

else
    nbond_types=0;
end

if nAngles>0

    angle_triplets=[[atom(Angle_index(:,1)).type]' [atom(Angle_index(:,2)).type]' [atom(Angle_index(:,3)).type]'];
    angle_info1=[];% angle_info2=[];
    for i = 1:size(angle_triplets,1)
        angle_info1{i,1} = sprintf('%s %s %s %.2f', angle_triplets{i,1}, angle_triplets{i,2}, angle_triplets{i,3}, Angle_index(i,4));
    end
    a1=angle_info1;[a1,idxa1]=unique(a1,'stable');
    nangle_types=size(a1,1);
    for i=1:size(angle_info1,1)
        [ind,angle_types(i)]=ismember(angle_info1(i),a1);
    end
    angle_types(angle_types>nangle_types)=angle_types(angle_types>nangle_types)-nangle_types;
    selected_Angle_index=Angle_index(idxa1,1:3);
    [H_row,H_col]=ind2sub(size(selected_Angle_index),find(ismember(selected_Angle_index,ind_H)));
    % Angle_index(H_row,1:4)

    angle_coeffs=[1:length(idxa1)]';
    angle_ka=Angle_index(idxa1,4);
    angle_deg=Angle_index(idxa1,4);
    angle_ka(:)=KANGLE;
    angle_ka(H_row)=KANGLEH;
    angle_deg(H_row)=angleH;
    angle_coeffs=[angle_coeffs angle_ka angle_deg];

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

% Start printing the lammps data file
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
fprintf(fid, '\n');
fprintf(fid, 'Masses \n');
fprintf(fid, '\n');

for i =1:natom_labels
    masses(i,:) = {i+prev_atom_types, num2str(Masses(i),precision), '#', Atom_labels{i}};
    fprintf(fid, '%i %-s %s %s\n', masses{i,:});
end


fprintf(fid, '\n');
%%

if size(ff,1)>0
    %%
    fprintf(fid, 'Pair Coeffs \n');
    fprintf(fid, '\n');
    for i =1:natom_labels
        paircoeffs(i,:) = {i, Epsilon(i), Sigma(i), '#', Atom_labels{i}};
        fprintf(fid, '%i %f %f %s %s\n', paircoeffs{i,:});
    end
else
    disp('Could not write pair coeffs!!')
end
fprintf(fid, '\n');
%
fprintf(fid, 'Atoms \n');
fprintf(fid, '\n');

for i = 1:nAtoms
    if sum(ismember(Atom_labels,[atom(i).type])) > 0
        Atom_label_ID(i,1)=find(ismember(Atom_labels,[atom(i).type])==1);
    end
    Atoms_data(i,:) = {i+prev_atom_index, [atom(i).molid]+prev_mol_index, Atom_label_ID(i,1), [atom(i).charge],[atom(i).x],[atom(i).y],[atom(i).z],'#',char(atom(i).type)};
    fprintf(fid, '\t%-i\t%-i\t%-i\t%8.5f\t%8.5f\t%8.5f\t%8.5f  %-s  %-s\n', Atoms_data{i,:});
end

%%
fprintf(fid, '\n');
fprintf(fid, '\n');



if nBonds>0

    %%
    fprintf(fid, 'Bond Coeffs \n');
    fprintf(fid, '\n');
    for i =1:size(bond_coeffs,1)
        temp_bond_string={bond_coeffs(i,:) strcat("# ",b1(i))};
        fprintf(fid, '\t%-i\t%-8.2f\t%-6.4f\t\t%s\n',temp_bond_string{:});
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    % Prints bond data
    fprintf(fid, 'Bonds \n');
    fprintf(fid, '\n');
    count_b = 1;
    while count_b <= nBonds
        Bond_order(count_b,:)={count_b+prev_bond_num, bond_types(count_b)+prev_bond_types, Bond_index(count_b,1)+prev_atom_index, Bond_index(count_b,2)+prev_atom_index,strcat("# ",[atom(Bond_index(count_b,1)).type],"-",[atom(Bond_index(count_b,2)).type]," ",num2str(Bond_index(count_b,3),3))};
        fprintf(fid, '\t%-i\t\t%-i\t\t%-i\t\t%-i\t\t\t%s\n', Bond_order{count_b,:});
        count_b = count_b + 1;
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');
end

if nAngles>0

    %%
    fprintf(fid, 'Angle Coeffs \n');
    fprintf(fid, '\n');
    for i =1:size(angle_coeffs,1)
        temp_angle_string={angle_coeffs(i,:) strcat("0      0       # ",a1(i))};
        fprintf(fid, '\t%-i\t\t%-8.2f\t%-8.2f\t\t%s\n',temp_angle_string{:});
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    % Prints angle data
    fprintf(fid, 'Angles \n');
    fprintf(fid, '\n');
    count_a = 1;
    while count_a <= nAngles % CHARMM st
        Angle_order(count_a,:)= {count_a+prev_angle_num, angle_types(count_a)+prev_angle_types, Angle_index(count_a,1)+prev_atom_index,Angle_index(count_a,2)+prev_atom_index,Angle_index(count_a,3)+prev_atom_index,strcat("# ",[atom(Angle_index(count_a,1)).type],"-",[atom(Angle_index(count_a,2)).type],"-",[atom(Angle_index(count_a,3)).type]," ",num2str(round(Angle_index(count_a,4),2)))};
        fprintf(fid, '\t%-i\t\t%-i\t\t%-i\t\t%-i\t\t%-i\t\t\t%s\n', Angle_order{count_a,:});
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
        fprintf(fid, '\t%-i\t%-i\t%-i\t%-i\t%-i\t\t%-i\n', Dihedral_order{count_d,:});
        count_d = count_d + 1;
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');
end

fclose(fid);

