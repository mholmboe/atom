%% write_minff_psf.m
% * This script creates and prints a .psf file
% * Works for minff.
% * The variables explicit_bonds|explicit_angles (1/0) on that are set
% * between lines ~50-120 for each specific forcefield, determines if the
% * bond and angle terms are added to the .itp file.
% *
% * In the examples below, the first cutoff (1.25) represents max bond
% * distance to any H. The second cutoff (2.45) represents the max bond
% * distance between any non-H atomtypes, like Si-O.
% * Additional commands governing the selection of bonds/angles can be
% * found on lines ~140-175, and 180-200 for angles
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_psf(atom,Box_dim,filename) % Basic input arguments
% # write_atom_psf(atom,Box_dim,filename,1.25,2.25,'minff','opc3')
%
function write_minff_psf(atom,Box_dim,filename,varargin)

format long
nAtoms=size(atom,2);

if regexp(filename,'.psf') ~= false
    filename = filename;
else
    filename = strcat(filename,'.psf');
end

if nargin > 3
    maxrshort=varargin{1}
    maxrlong=varargin{2}
else
    maxrshort=1.25;
    maxrlong=2.45;
end

if nargin>5
    ffname=varargin{3}
    if nargin>6
        watermodel=varargin{4}
    else
        disp('Unknown watermodel, will try OPC3')
        watermodel='OPC3'
    end
end
disp('Forcefield not stated, will make some assumptions then...')
ffname='minff';
watermodel='OPC3'; % SPC/E, depreceated
Atom_labels=unique([atom.type]);
atom = mass_atom(atom);
Masses=[];
for i=1:size(Atom_labels,2)
    ind=find(strcmpi([atom.type],Atom_labels(i)));
    Masses(i)=[atom(ind(i)).mass];
end

if ~isfield(atom,'charge')
    atom=charge_minff_atom(atom,Box_dim,{'Al' 'Alt' 'Ale' 'Tio' 'Feo' 'Fet' 'Fee' 'Fe2' 'Fe2e' 'Fe3e' 'Na' 'K' 'Cs' 'Mgo' 'Mgh' 'Mge' 'Cao' 'Cah' 'Sit' 'Si' 'Sio' 'Site' 'Lio' 'H'},[1.782 1.782 1.985 2.48 1.14 1.14 1.14 0.7 0.86666 1.45 1 1 1 1.562 1.74 1.635 1.66 1.52 1.884 1.884 1.884 2.413 0.86 0.4]);
end
Total_charge=sum([atom.charge])
round2dec(Total_charge,5)
%         pause
nrexcl=1; % See the gromacs manual
explicit_bonds = 0
explicit_angles = 1


%% Find atomtype specific indexes

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

%
file_title = 'An awesome psf file'; % Header in output file
molecule_name = char([atom(1).resname]); % molecule name
Atom_label = unique([atom.type]);

fid = fopen(filename, 'wt');

fprintf(fid, '%-s\n','PSF');
fprintf(fid, '\n');
fprintf(fid, '%s\n','       2 !NTITLE');
fprintf(fid, '%s\n',' REMARKS MATLAB-generated PSF structure file');
fprintf(fid, '%s\n',' REMARKS coded by MHolmboe (michael.holmboe@umu.se)');
fprintf(fid, '\n');
fprintf(fid, '%8i %s\n',nAtoms,'!NATOM');

XYZ_labels=[atom.type];
Atom_label=unique([atom.type]);
atomID = 1:size([atom.type],2);
molID=zeros(1,size([atom.type],2));
Atom_label_ID=zeros(size([atom.type],2),1);

for i = 1:length(Atom_label)
    Atom_label_ID(ismember([atom.type],Atom_label(i)))=i;
end

assignin('base','atom_psf',atom);

ResNum=[atom.molid];
SegName='MOL';
ResName=[atom.resname];
for i = 1:size([atom.type],2)
    Atoms_data(i,:) = {atomID(1,i),char(SegName),[atom(i).molid],char([atom(i).resname]),char([atom(i).element]),char([atom(i).type]),[atom(i).charge],[atom(i).mass],0};
    fprintf(fid,'%8i %4s %-5i%-5s%-5s%-5s%10.6f%14.4f%12i\n', Atoms_data{i,:});
end

fprintf(fid, '\n');
fprintf(fid, '%8i %s\n',nBonds,'!NBOND: bonds');

bond_matrix=Bond_index(:,1:2);
while mod(2*size(bond_matrix,1),8) ~= 0
    bond_matrix=[bond_matrix;0 0];
end
bond_temp=zeros(2*size(bond_matrix,1),1);
bond_temp(1:2:end,1)=bond_matrix(:,1);
bond_temp(2:2:end,1)=bond_matrix(:,2);
bond_list=[reshape(bond_temp',8,[])]';

count_b = 1;
bondtype=1;
while count_b <= size(bond_list,1)%length(bond_list)
    if bond_list(count_b,1)~=0;fprintf(fid, '%8i%8i',bond_list(count_b,1:2));end
    if bond_list(count_b,3)~=0;fprintf(fid, '%8i%8i',bond_list(count_b,3:4));end
    if bond_list(count_b,5)~=0;fprintf(fid, '%8i%8i',bond_list(count_b,5:6));end
    if bond_list(count_b,7)~=0;fprintf(fid, '%8i%8i',bond_list(count_b,7:8));end
    fprintf(fid, '\n');
    count_b = count_b + 1;
end

fprintf(fid, '\n');

fprintf(fid, '%8i %s\n',nAngles,'!NTHETA: angles');
angle_matrix=Angle_index(:,1:3);
while mod(3*size(angle_matrix,1),9) ~= 0
    angle_matrix=[angle_matrix;0 0 0];
end
angle_temp=zeros(3*size(angle_matrix,1),1);
angle_temp(1:3:end)=angle_matrix(:,1);
angle_temp(2:3:end)=angle_matrix(:,2);
angle_temp(3:3:end)=angle_matrix(:,3);
angle_list=[reshape(angle_temp',9,[])]';

count_a = 1;
angletype=1;
while count_a <= size(angle_list,1)
    if angle_list(count_a,1)~=0;fprintf(fid, '%8i%8i%8i',angle_list(count_a,1:3));end
    if angle_list(count_a,4)~=0;fprintf(fid, '%8i%8i%8i',angle_list(count_a,4:6));end
    if angle_list(count_a,7)~=0;fprintf(fid, '%8i%8i%8i',angle_list(count_a,7:9));end
    fprintf(fid, '\n');
    count_a = count_a + 1;
end


%
%
%        0 !NPHI: dihedrals
%
%
%        0 !NIMPHI: impropers
%
%
%        0 !NDON: donors
%
%
%        0 !NACC: acceptors
%
%
%        0 !NNB
%
%
%        1       0 !NGRP
%        0       0       0

fprintf(fid, '\n');

fprintf(fid, '\n');
fprintf(fid, '%8i %s\n',0,'!NPHI: dihedrals');
fprintf(fid, '\n');

fprintf(fid, '\n');
fprintf(fid, '%8i %s\n',0,'!NIMPHI: impropers');
fprintf(fid, '\n');

donor_ind=sort(unique([find(strncmp([atom.type],'O',1)) find(strncmp([atom.type],'N',1))]));

fprintf(fid, '\n');
fprintf(fid, '%8i %s\n',0,'!NDON: donors');
fprintf(fid, '\n');
acceptor_ind=sort(unique([find(strncmp([atom.type],'O',1)) find(strncmp([atom.type],'N',1))]));

fprintf(fid, '\n');
fprintf(fid, '%8i %s\n',0,'!NACC: acceptors');
fprintf(fid, '\n');

fprintf(fid, '\n');
fprintf(fid, '%8i %s\n',0,'!NNB');
fprintf(fid, '\n');

fprintf(fid, '\n');
fprintf(fid, '%8i %8i %s\n',1,0,'!NGRP');
fprintf(fid, '%8i %8i %8i\n',0,0,0);
fprintf(fid, '\n');

fclose(fid);

assignin('caller','psf_atom',atom);

