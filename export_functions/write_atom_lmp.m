%% write_atom_lmp.m
% * This script creates and prints a lammps data file (.data). Works best for
% Clayff systems, defining the ffname 'clayff' and watermodel 'spce'.
% Nevertheless, this new version should be able to handle bonds|angles|dihedrals
%
%
%% Version
% 3.00
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
ind_OW=strncmpi([atom.type],'Ow',2);[atom(ind_OW).type]=deal({'Ow'});
ind_HW=strncmpi([atom.type],'Hw',2);[atom(ind_HW).type]=deal({'Hw'});
Atom_labels=sort(unique([atom.type]));
natom_labels=size(Atom_labels,2);

if regexp(filename,'.data') ~= false
    filename = filename;
else
    filename = strcat(filename,'.data');
end

if nargin > 3
    short_r=varargin{1};
    long_r=varargin{2};
else
    short_r=1.25;
    long_r=2.45; % long=short since clayff
end

% if nargin>5
%     ffname=varargin(3);
%     if nargin>6
%         watermodel=varargin(4)
%     else
%         disp('Unknown watermodel, will try SPC/E')
%         watermodel='SPC/E';
%     end
%     if strcmpi(ffname,'clayff_2004')
%         clayff_2004_param(sort(unique([atom.type])),watermodel);
%         if ~isfield(atom,'charge')
%             atom = charge_atom(atom,Box_dim,'clayff_2004',watermodel,'adjust');
%         end
%         Total_charge=sum([atom.charge])
%         round2dec(Total_charge,5)
%     elseif strncmpi(ffname,'clayff',5)
%         clayff_param(sort(unique([atom.type])),watermodel);
%         atom = charge_atom(atom,Box_dim,'clayff',watermodel,'more');
%         Total_charge
%         round2dec(Total_charge,5)
%     end
% else
%     %     disp('Unknown forcefield, will try clayff')
%     %     clayff_param(sort(unique([atom.type])),'spc/e');
%     %     atom = charge_atom(atom,Box_dim,'clayff','spc/e');
%     %     Total_charge
%     disp('Forcefield not stated, will make some assumptions then...')
%     pause(2)
%     ffname='clayff_2004'
%     watermodel='SPC/E'
%     clayff_2004_param(sort(unique([atom.type])),watermodel);
%     if ~isfield(atom,'charge')
%         atom = charge_atom(atom,Box_dim,'clayff_2004',watermodel,'adjust');
%     end
%     Total_charge=sum([atom.charge])
%     round2dec(Total_charge,5)
% end

if nargin>5
    ffname=varargin{3}
    if nargin>6
        watermodel=varargin{4}
    else
        disp('Unknown watermodel, will try SPC/E')
        watermodel='SPC/E'
    end

    if strcmpi(ffname,'minff')
        atom = mass_atom(atom);
        Total_charge=sum([atom.charge])
        round2dec(Total_charge,5)
        %         pause
        nrexcl=1; % See the gromacs manual
        explicit_bonds = 0
        explicit_angles = 1
        watermodel='OPC3'; % SPC/E, depreceated
    elseif strcmpi(ffname,'clayff')
        clayff_param(sort(unique([atom.type])),watermodel);
        if ~isfield(atom,'charge')
            atom = charge_atom(atom,Box_dim,'clayff',watermodel,'adjust');
        end
        Total_charge=sum([atom.charge])
        round2dec(Total_charge,5)
        %         pause
        nrexcl=1; % See the gromacs manual
        explicit_bonds = 0
        explicit_angles = 1
    elseif strncmpi(ffname,'clayff_2004',5)
        clayff_2004_param(sort(unique([atom.type])),watermodel);
        if ~isfield(atom,'charge')
            atom = charge_atom(atom,Box_dim,'clayff_2004',watermodel,'adjust');
        end
        Total_charge=sum([atom.charge])
        round2dec(Total_charge,5)
        %         pause
        nrexcl=1; % See the gromacs manual
        explicit_bonds = 0
        explicit_angles = 1
    elseif strcmpi(ffname,'interface')
        interface_param(sort(unique([atom.type])),watermodel);
        if ~isfield(atom,'charge')
            atom = charge_atom(atom,Box_dim,'interface',watermodel,'adjust');
        end
        Total_charge=sum([atom.charge])
        nrexcl=3; % See the gromacs manual
        explicit_bonds = 1;
        explicit_angles = 1;
    elseif strcmpi(ffname,'interface15')
        atom = mass_atom(atom);
        nrexcl=3; % See the gromacs manual
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
        explicit_bonds = 0 % 0 currently does not work, because no default bond types
        explicit_angles = 0% 0 currently does not work, because no default angle types
    elseif strcmpi(ffname,'interface_car')
        % Experimental!!!
        atom = mass_atom(atom);
        nrexcl=3; % See the gromacs manual
        explicit_bonds = 1
        explicit_angles = 1
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
    ffname='minff';
    watermodel='OPC3'; % SPC/E, depreceated
    % minff_param(sort(unique([atom.type])),'OPC3');
    atom = mass_atom(atom);
    if ~isfield(atom,'charge')
        atom=charge_minff_atom(atom,Box_dim,{'Al' 'Alt' 'Ale' 'Tio' 'Feo' 'Fet' 'Fee' 'Fe3e' 'Fe2' 'Fe2e' 'Na' 'K' 'Cs' 'Mgo' 'Mgh' 'Mge' 'Cao' 'Cah' 'Sit' 'Si' 'Sio' 'Site' 'Lio' 'H'},[1.782 1.782 1.985 2.48 1.5 1.5 1.75 1.75 1.184 1.32 1 1 1 1.562 1.74 1.635 1.66 1.52 1.884 1.884 1.884 2.413 0.86 0.4]);
    end
    Total_charge=sum([atom.charge])
    round2dec(Total_charge,5)
    %         pause
    nrexcl=1; % See the gromacs manual
    explicit_bonds = 0
    explicit_angles = 1
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

    %% To reduce the number of bond types
    bond_pairs=sort([[atom(Bond_index(:,1)).type]' [atom(Bond_index(:,2)).type]']);

    bond_info1=[];
    for i = 1:size(bond_pairs,1)
        
        bond_info1{i,1} = sprintf('%s %s', bond_pairs{i,1}, bond_pairs{i,2});
    end
    b1=bond_info1;[b1,idxb1]=unique(b1,'stable');
    nbond_types=size(b1,1);
    for i=1:size(bond_info1,1)
        [ind,bond_types(i)]=ismember(bond_info1(i),[b1]); %;b2]);
    end
    bond_types(bond_types>nbond_types)=bond_types(bond_types>nbond_types)-nbond_types;

    bond_coeffs=[1:length(idxb1)]';
    bond_dist=Bond_index(idxb1,3);
    bond_kb=Bond_index(idxb1,3);
    bond_kb(bond_dist>1.25)=kbM;
    bond_kb(bond_dist<1.25)=kbH;

    bond_dist(bond_dist<1.25)=bHdist;
    bond_coeffs=[bond_coeffs bond_kb bond_dist];

    % %% To maximize the number of bond types
    % bond_pairs=[[atom(Bond_index(:,1)).type]' [atom(Bond_index(:,2)).type]'];
    % bond_info1=[];%bond_info2=[];
    % for i = 1:size(bond_pairs,1)
    %     bond_info1{i,1} = sprintf('%s %s %.3f', bond_pairs{i,1}, bond_pairs{i,2}, Bond_index(i,3));
    %     % bond_info2{i,1} = sprintf('%s %s %.3f', bond_pairs{i,2}, bond_pairs{i,1}, Bond_index(i,3));
    % end
    % b1=bond_info1;[b1,idxb1]=unique(b1,'stable');
    % % b2=bond_info1;b2=unique(b2,'stable');
    % nbond_types=size(b1,1);
    % for i=1:size(bond_info1,1)
    %     [ind,bond_types(i)]=ismember(bond_info1(i),[b1]); %;b2]);
    %     % bond_coeffs
    % end
    % bond_types(bond_types>nbond_types)=bond_types(bond_types>nbond_types)-nbond_types;
    % 
    % bond_coeffs=[1:length(idxb1)]';
    % bond_dist=Bond_index(idxb1,3);
    % bond_kb=Bond_index(idxb1,3);
    % bond_kb(bond_dist>1.25)=kbM;
    % bond_kb(bond_dist<1.25)=kbH;
    % 
    % bond_dist(bond_dist<1.25)=bHdist;
    % bond_coeffs=[bond_coeffs bond_kb bond_dist];

else
    nbond_types=0;
end

if nAngles>0

     %% To reduce the number of angle types
    angle_triplets=[[atom(Angle_index(:,1)).type]' [atom(Angle_index(:,2)).type]' [atom(Angle_index(:,3)).type]';...
        [atom(Angle_index(:,3)).type]' [atom(Angle_index(:,2)).type]' [atom(Angle_index(:,1)).type]'];
    Angle_index2=[Angle_index(:,1:4);Angle_index(:,[3 2 1 4])];

    
    angle_endtypes=string(angle_triplets(:,[1,3]));
    angle_endtypes(:, [1, 2]) = cellstr(angle_endtypes);
    angle_triplets=[angle_endtypes(:,1) [atom(Angle_index2(:,2)).type]' angle_endtypes(:,2)];

    angle_info1=[];% angle_info2=[];
    for i = 1:size(angle_triplets,1)
        angle_info1{i,1} = sprintf('%s %s %s', angle_triplets{i,1}, angle_triplets{i,2}, angle_triplets{i,3});
    end
    a1=angle_info1;[a1,idxa1]=unique(a1,'stable');
    nangle_types=size(a1,1);
    for i=1:size(angle_info1,1)
        [ind,angle_types(i)]=ismember(angle_info1(i),a1);
    end
    angle_types(angle_types>nangle_types)=angle_types(angle_types>nangle_types)-nangle_types;
    selected_Angle_index=Angle_index2(idxa1,1:3);
    [H_row,H_col]=ind2sub(size(selected_Angle_index),find(ismember(selected_Angle_index,ind_H)));
    [Hw_row,Hw_col]=ind2sub(size(selected_Angle_index),find(ismember(selected_Angle_index,ind_Hw)));
    % Angle_index(H_row,1:4)

    angle_coeffs=[1:length(idxa1)]';
    for i=1:size(idxa1,1)
        angle_deg(i,1)=mean(Angle_index2(ismember(angle_info1,angle_info1(idxa1(i),:)),4));
    end
    
    angle_ka=angle_deg;
    angle_ka(:)=KANGLE;

    angle_ka(H_row)=KANGLEH;
    angle_deg(H_row)=angleH;

    angle_ka(Hw_row)=KANGLE_WAT;
    angle_deg(Hw_row)=ANGLE_WAT;

    angle_coeffs=[angle_coeffs angle_ka angle_deg];

    % %% To maximize the number of angles
    % 
    % angle_triplets=[[atom(Angle_index(:,1)).type]' [atom(Angle_index(:,2)).type]' [atom(Angle_index(:,3)).type]'];
    % angle_info1=[];% angle_info2=[];
    % for i = 1:size(angle_triplets,1)
    %     angle_info1{i,1} = sprintf('%s %s %s %.2f', angle_triplets{i,1}, angle_triplets{i,2}, angle_triplets{i,3}, Angle_index(i,4));
    % end
    % a1=angle_info1;[a1,idxa1]=unique(a1,'stable');
    % nangle_types=size(a1,1);
    % for i=1:size(angle_info1,1)
    %     [ind,angle_types(i)]=ismember(angle_info1(i),a1);
    % end
    % angle_types(angle_types>nangle_types)=angle_types(angle_types>nangle_types)-nangle_types;
    % selected_Angle_index=Angle_index(idxa1,1:3);
    % [H_row,H_col]=ind2sub(size(selected_Angle_index),find(ismember(selected_Angle_index,ind_H)));
    % % Angle_index(H_row,1:4)
    % 
    % angle_coeffs=[1:length(idxa1)]';
    % angle_ka=Angle_index(idxa1,4);
    % angle_deg=Angle_index(idxa1,4);
    % angle_ka(:)=KANGLE;
    % angle_ka(H_row)=KANGLEH;
    % angle_deg(H_row)=angleH;
    % angle_coeffs=[angle_coeffs angle_ka angle_deg];

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

% Start printing the lammps .data file
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

if isfield(atom,'mass')
    for i =1:natom_labels
        masses(i,:) = {i+prev_atom_types, num2str([atom(i).mass],precision), '#', Atom_labels{i}};
        fprintf(fid, '%i %-s %s %s\n', masses{i,:});
    end
elseif exist('Masses','var')
    for i =1:natom_labels
        masses(i,:) = {i+prev_atom_types, num2str(Masses(i),precision), '#', Atom_labels{i}};
        fprintf(fid, '%i %-s %s %s\n', masses{i,:});
    end
end

fprintf(fid, '\n');
%%
%%
fprintf(fid, 'Pair Coeffs \n');
fprintf(fid, '\n');

% for i =1:natom_labels
%     paircoeffs(i,:) = {i, Epsilon(i)*kcalmol_to_eV, Sigma(i), '#', Atom_labels{i}};
%     fprintf(fid, '%i %E %f %s %s\n', paircoeffs{i,:});
% end
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

