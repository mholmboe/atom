%% write_atom_psf.m
% * This function writes an .psf file from the atom struct
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_psf(atom,Box_dim,filename_out) % Not recommended, Please specify which cutoff's and forcefield to use
% # write_atom_psf(atom,Box_dim,filename_out,1.25,1.25,'clayff','spce')
%
function write_atom_psf(atom,Box_dim,filename,varargin)

nAtoms=size(atom,2);

if regexp(filename,'.psf') ~= false
    filename = filename;
else
    filename = strcat(filename,'.psf');
end

Atom_label = unique([atom.type]);
if nargin > 3
    short_r=cell2mat(varargin(1));
    long_r=cell2mat(varargin(2));
else
    short_r=1.25;
    long_r=2.25;
end

if nargin>5
    ffname=varargin(3);
    if nargin>6
        watermodel=varargin(4);
    else
        disp('Unknown watermodel, will try SPC/E')
        watermodel='SPC/E';
    end
    if strncmpi(ffname,'clayff',5)
        clayff_param(sort(unique([atom.type])),watermodel);
        if ~isfield(atom,'charge')
            atom = charge_atom(atom,Box_dim,'clayff',watermodel,'adjust');
        end
        Total_charge = check_clayff_charge(atom)
    elseif strncmpi(ffname,'clayff_2004',5)
        clayff_2004_param(sort(unique([atom.type])),watermodel);
        if ~isfield(atom,'charge')
            atom = charge_atom(atom,Box_dim,'clayff_2004',watermodel,'adjust');
        end
        Total_charge=sum([atom.charge])
    elseif strcmpi(ffname,'interface')
        interface_param(sort(unique([atom.type])),watermodel);
        if ~isfield(atom,'charge')
            atom = charge_atom(atom,Box_dim,'interface','tip3p');
        end
        Total_charge = check_interface_charge(atom)
    elseif strcmpi(ffname,'interface15')
        if ~isfield(atom,'charge')
            atom = charge_atom(atom,Box_dim,'interface15','tip3p');
        end
        atom = check_interface15_charge(atom,'CLAY_MINERALS');
        atom = mass_atom(atom);
    elseif strcmpi(ffname,'interface_car')
        % Experimental!!!
        atom = mass_atom(atom);
        nrexcl=2; % See the gromacs manual
        explicit_bonds = 0;
        explicit_angles = 0;
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
    %         clayff_param(sort(unique([atom.type])),watermodel);
    %         Total_charge = check_clayff_charge(atom)
end

lx=Box_dim(1);ly=Box_dim(2);lz=Box_dim(3);
if length(Box_dim)>3
    xy=Box_dim(6);xz=Box_dim(8);yz=Box_dim(9);
else
    xy=0;xz=0;yz=0;
end

fid = fopen(filename, 'wt');

fprintf(fid, '%-s\r\n','PSF');
fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','       2 !NTITLE');
fprintf(fid, '%s\r\n',' REMARKS MATLAB-generated PSF structure file');
fprintf(fid, '%s\r\n',' REMARKS coded by MHolmboe (michael.holmboe@umu.se)');
fprintf(fid, '\r\n');
fprintf(fid, '%8i %s\r\n',nAtoms,'!NATOM');

XYZ_labels=[atom.type];
Atom_label=unique([atom.type]);
bond_angle_atom(atom,Box_dim,short_r,long_r);
atomID = 1:size([atom.type],2);
molID=zeros(1,size([atom.type],2));
Atom_label_ID=zeros(size([atom.type],2),1);

% if exist('ffname','var')
%     atom = charge_atom(atom,Box_dim,ffname,watermodel);
% end

for i = 1:length(Atom_label)
    Atom_label_ID(ismember([atom.type],Atom_label(i)))=i;
end

ResNum=[atom.molid];
SegName=['SURF'];
ResName=[atom.resname];
for i = 1:size([atom.type],2)
    
    if exist('Masses','var')
        Atom_label_ID(i,1);Charge(Atom_label_ID(i,1));
        %                 atomID,     segname, residueID,  resname,       atomname,                      atomtype,                      charge,                     mass,        and an unused 0
        Atoms_data(i,:) = [atomID(1,i),SegName,[atom(i).molid],[atom(i).resname],char([atom(i).type]),char([atom(i).type]),[atom(i).charge],Masses(Atom_label_ID(i,1)),0];
    else
        i;
        [atom(i).mass];
        [atom(i).charge];
        Atoms_data(i,:) = [atomID(1,i),SegName,[atom(i).molid],[atom(i).resname],char([atom(i).element]),char([atom(i).type]),[atom(i).charge],[atom(i).mass],0];
    end
    
    fprintf(fid, '%8i %4s %-5i%-5s%-5s%-5s%10.6f%14.4f%12i\r\n', Atoms_data{i,:});
end

fprintf(fid, '\r\n');
fprintf(fid, '%8i %s\r\n',nBonds,'!NBOND: bonds');

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
while count_b <= length(bond_list)
    if bond_list(count_b,1)~=0;fprintf(fid, '%8i%8i',bond_list(count_b,1:2));end
    if bond_list(count_b,3)~=0;fprintf(fid, '%8i%8i',bond_list(count_b,3:4));end
    if bond_list(count_b,5)~=0;fprintf(fid, '%8i%8i',bond_list(count_b,5:6));end
    if bond_list(count_b,7)~=0;fprintf(fid, '%8i%8i',bond_list(count_b,7:8));end
    fprintf(fid, '\r\n');
    count_b = count_b + 1;
end

fprintf(fid, '\r\n');

fprintf(fid, '%8i %s\r\n',nAngles,'!NTHETA: angles');
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
    fprintf(fid, '\r\n');
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

fprintf(fid, '\r\n');

fprintf(fid, '\r\n');
fprintf(fid, '%8i %s\r\n',0,'!NPHI: dihedrals');
fprintf(fid, '\r\n');

fprintf(fid, '\r\n');
fprintf(fid, '%8i %s\r\n',0,'!NIMPHI: impropers');
fprintf(fid, '\r\n');

donor_ind=sort(unique([find(strncmp([atom.type],'O',1)) find(strncmp([atom.type],'N',1))]));

fprintf(fid, '\r\n');
fprintf(fid, '%8i %s\r\n',0,'!NDON: donors');
fprintf(fid, '\r\n');
acceptor_ind=sort(unique([find(strncmp([atom.type],'O',1)) find(strncmp([atom.type],'N',1))]));

fprintf(fid, '\r\n');
fprintf(fid, '%8i %s\r\n',0,'!NACC: acceptors');
fprintf(fid, '\r\n');

fprintf(fid, '\r\n');
fprintf(fid, '%8i %s\r\n',0,'!NNB');
fprintf(fid, '\r\n');

fprintf(fid, '\r\n');
fprintf(fid, '%8i %8i %s\r\n',1,0,'!NGRP');
fprintf(fid, '%8i %8i %8i\r\n',0,0,0);
fprintf(fid, '\r\n');

fclose(fid);

assignin('caller','itp_atom',atom);

