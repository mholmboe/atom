%% write_atom_psf.m
% * This function writes an .psf file from the atom struct
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_psf(atom,Box_dim,filename_out) % Basic input arguments
% # write_atom_psf(atom,Box_dim,filename_out,1.25,1.25,'clayff','spce') % Specifying the rmaxshort and rmaxlong cutoff's and forcefield to use
%
function write_atom_psf(atom,Box_dim,filename,varargin)

format long
nAtoms=size(atom,2);

if regexp(filename,'.psf') ~= false
    filename = filename;
else
    filename = strcat(filename,'.psf');
end

Atom_label = unique([atom.type]);

if nargin > 3
    maxrshort=cell2mat(varargin(1))
    maxrlong=cell2mat(varargin(2))
else
    maxrshort=1.25;
    maxrlong=1.25; % long=short since clayff
end

if nargin>5
    ffname=varargin(3)
    if nargin>6
        watermodel=varargin(4)
    else
        disp('Unknown watermodel, will try SPC/E')
        watermodel='SPC/E'
    end
    
    if strcmpi(ffname,'clayff')
        clayff_param(sort(unique([atom.type])),watermodel);
        atom=mass_atom(element_atom(atom));
        if ~isfield(atom,'charge')
            atom = charge_atom(atom,Box_dim,'clayff',watermodel,'adjust');
        end
        Total_charge=sum([atom.charge])
        round(Total_charge,5)
        %         pause
        nrexcl=1; % See the gromacs manual
        explicit_bonds = 0;
        explicit_angles = 1;
    elseif strncmpi(ffname,'clayff_2004',5)
        clayff_2004_param(sort(unique([atom.type])),watermodel);
        if ~isfield(atom,'charge')
            atom = charge_atom(atom,Box_dim,'clayff_2004',watermodel,'adjust');
        end
        Total_charge=sum([atom.charge])
        round(Total_charge,5)
        %         pause
        nrexcl=1; % See the gromacs manual
        explicit_bonds = 0;
        explicit_angles = 1;
    elseif strcmpi(ffname,'interface')
        interface_param(sort(unique([atom.type])),watermodel);
        if ~isfield(atom,'charge')
            atom = charge_atom(atom,Box_dim,'interface',watermodel,'adjust');
        end
        Total_charge=sum([atom.charge])
        nrexcl=2; % See the gromacs manual
        explicit_bonds = 1;
        explicit_angles = 1;
        
    elseif strcmpi(ffname,'interface15')
        atom = mass_atom(atom);
        nrexcl=2; % See the gromacs manual
        %         interface15_param(sort(unique([atom.type])),watermodel);
        %         atom = charge_atom(atom,Box_dim,'interface15',watermodel,'adjust');
        if nargin > 7
            model_database=varargin{5}
        else
            model_database='CLAY_MINERALS'
        end
        atom = check_interface15_charge(atom,model_database);
        Total_charge
        nrexcl=2; % See the gromacs manual
        explicit_bonds = 0; % 0 currently does not work, because no default bond types
        explicit_angles = 0;% 0 currently does not work, because no default angle types
    elseif strcmpi(ffname,'interface_car')
        % Experimental!!!
        atom = mass_atom(atom);
        nrexcl=2; % See the gromacs manual
        explicit_bonds = 1;
        explicit_angles = 1;
        %     elseif strcmpi(ffname,'oplsaa_go');
        %         % This is not for you...
        %         oplsaa_go_param(sort(unique([atom.type])),watermodel);
        %         atom = charge_opls_go_atom(atom,Box_dim,{'H' 'Oe' 'Oh'},[0.418 -0.4 -0.683])
        %         Total_charge
        %         nrexcl=3; % See the gromacs manual
        %         explicit_bonds = 0;
        %         explicit_angles = 0;
    end
else
    disp('Forcefield not stated, will make some assumptions then...')
    pause(2)
    ffname='clayff_2004'
    watermodel='SPC/E'
    pause(2)
    atom = mass_atom(atom);
    element=element_atom(atom);
    [atom.element]=element.type;
    if ~isfield(atom,'charge')
        atom = charge_atom(atom,Box_dim,ffname,watermodel);
    end
    nrexcl=1; % See the gromacs manual
    explicit_bonds = 0;
    explicit_angles = 0;
end

%% Find atomtype specific indexes

ind_Hneighbours = find(~cellfun(@isempty,regexpi([atom.type],'h')));
ind_H=find(strncmpi([atom.type],{'H'},1));
ind_O=find(strncmpi([atom.type],{'O'},1));
ind_Osih=find(strncmpi([atom.type],{'Osih'},4));
ind_Oh=intersect(ind_O,ind_Hneighbours);
ind_Al=find(strncmpi([atom.type],'Al',2));
ind_Mgo=find(strncmpi([atom.type],{'Mgo'},3));
ind_Si=find(strncmpi([atom.type],{'Si'},2));
ind_Oct=sort([ind_Mgo]);

atom = bond_angle_atom(atom,Box_dim,maxrshort,maxrlong);
% if strncmpi(ffname,'clayff',5)
%     %     %% To only keep bonds to atoms also bonded to H's, uncomment the next four lines
%     %     disp('Keeping only bonds with H')
%     %     [h_row,h_col]=ind2sub(size(Bond_index),find(ismember(Bond_index,ind_Hneighbours)));
%     %     Bond_index=Bond_index(h_row,:);
%     %     nBonds=size(Bond_index,1);
%
%% To only keep bonds to H's, uncomment the next three lines
[H_row,H_col]=ind2sub(size(Bond_index),find(ismember(Bond_index,ind_H)));
Bond_index=Bond_index(H_row,:);
nBonds=size(Bond_index,1);
%
%     %     %% To only keep bonds between Osih - H, uncomment the next four lines
%     %     disp('Keeping only bonds with H')
%     %     [h_row,h_col]=ind2sub(size(Bond_index),find(ismember(Bond_index,ind_Osih)));
%     %     Bond_index=Bond_index(h_row,:);
%     %     nBonds=size(Bond_index,1);
%
%     %     %% To remove bonds with 'Al'
%     %     [Al_row,Al_col]=ind2sub(size(Bond_index),find(ismember(Bond_index,ind_Al)));
%     %     Bond_index(Al_row,:)=[];
%     %     nBonds=size(Bond_index,1);
%
%     %     %% To remove bonds with 'Si'
%     %     [Si_row,Si_col]=ind2sub(size(Bond_index),find(ismember(Bond_index(:,2),ind_Si)));
%     %     Bond_index(Si_row,:)=[];
%     %     nBonds=size(Bond_index,1);
%
%     %    %% To remove bonds larger than certain rmin, uncomment next two lines
%     %     rm_ind=find(Bond_index(:,3)>1.25);
%     %     Bond_index(rm_ind,:)=[];
%     %     nBonds=size(Bond_index,1);
%
% end

[Y,I]=sort(Bond_index(:,1));
Bond_index=Bond_index(I,:);
Bond_index = unique(Bond_index,'rows','stable');

% if strncmpi(ffname,'clayff',5)
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
file_title = 'Gromacs awesome itp file'; % Header in output file
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
% bond_angle_atom(atom,Box_dim,short_r,long_r);
atomID = 1:size([atom.type],2);
molID=zeros(1,size([atom.type],2));
Atom_label_ID=zeros(size([atom.type],2),1);

% if exist('ffname','var')
%     atom = charge_atom(atom,Box_dim,ffname,watermodel);
% end

for i = 1:length(Atom_label)
    Atom_label_ID(ismember([atom.type],Atom_label(i)))=i;
end

assignin('base','atom_psf',atom);

ResNum=[atom.molid];
SegName='MOL';
ResName=[atom.resname];
for i = 1:size([atom.type],2)
    
    if isfield(atom,'mass')
        i;
        atomID(1,i)
        SegName
        [atom(i).molid]
        [atom(i).resname]
        char([atom(i).element])
        char([atom(i).type])
        [atom(i).charge]
        [atom(i).mass]
        Atoms_data(i,:) = {atomID(1,i),char(SegName),[atom(i).molid],char([atom(i).resname]),char([atom(i).element]),char([atom(i).type]),[atom(i).charge],[atom(i).mass],0};
    elseif exist('Masses','var')
        Masses
        Atom_label_ID(i,1);Charge(Atom_label_ID(i,1));
        %                 atomID,     segname, residueID,  resname,       atomname,                      atomtype,                      charge,                     mass,        and an unused 0
        Atoms_data(i,:) = {atomID(1,i),char(SegName),[atom(i).molid],char([atom(i).resname]),char([atom(i).type]),char([atom(i).type]),[atom(i).charge],Masses(Atom_label_ID(i,1)),0};
    end
    Atoms_data(i,:)
    Atoms_data{i,:}
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

assignin('caller','itp_atom',atom);

