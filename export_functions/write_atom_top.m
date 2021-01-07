%% write_atom_top.m
% * This script creates and prints a gromacs .top file
% * Works best for clayff or interface ff with spc, spce or tip3p
% * The variables explicit_bonds|explicit_angles (1/0) on that are set
% * between lines ~50-120 for each specific forcefield, determines if the
% * bond and angle terms are added to the .itp file.
% *
% * In the examples below, the first cutoff (1.25) represents max bond
% * distance to any H. The second cutoff (2.25) represents the max bond
% * distance between any non-H atomtypes, like Si-O.
% * Additional commands governing the selection of bonds/angles can be
% * found on lines ~140-170, and 175-190 for angles
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_itp(atom,Box_dim,filename) % Basic input arguments
% # write_atom_itp(atom,Box_dim,filename,1.25,1.25) % Default forcefield is clayff_2004
% # write_atom_itp(atom,Box_dim,filename,1.25,2.25,'clayff','spc/e')
% # write_atom_itp(atom,Box_dim,filename,1.25,2.25,'interface','tip3p')

function write_atom_top(atom,Box_dim,filename,varargin)

format long
nAtoms=size(atom,2);

if regexp(filename,'.itp') ~= false
    filename = filename;
else
    filename = strcat(filename,'.itp');
end

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
    %% To only keep bonds to H's, uncomment the next four lines
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

fid = fopen(filename, 'wt'); % open a text file that we can write into

fprintf(fid, '%s % s\r\n',';',file_title);
fprintf(fid, '%s % s\r\n',';','File written by MHolmboe (michael.holmboe@umu.se)');
fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','[ moleculetype ]');
fprintf(fid, '%s % s\r\n',';','molname   nrexcl');
% fprintf(fid, '%s       %d\r\n',strrep(molecule_name(1:3),'.itp',''),nrexcl);
fprintf(fid, '%s       %d\r\n',molecule_name(1:3),nrexcl);
fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','[ atoms ]');
fprintf(fid, '%s\r\n','; id   attype  resnr resname  atname   cgnr  charge      mass');

Atom_label_ID=ones(size(atom,2),1);
for i = 1:nAtoms
    if sum(ismember(Atom_label,[atom(i).type])) > 0
        Atom_label_ID(i,1)=find(ismember(Atom_label,[atom(i).type])==1);
    end
    if isfield(atom,'mass')
        Atoms_data(i,:) = {i, char([atom(i).fftype]),[atom(i).molid],molecule_name(1:3),char([atom(i).type]),i, round([atom(i).charge],6),[atom(i).mass]};
    else exist('Masses','var');
        Atoms_data(i,:) = {i, char([atom(i).type]),[atom(i).molid],molecule_name(1:3),char([atom(i).type]),i, round([atom(i).charge],6), Masses(Atom_label_ID(i,1))};
    end
    fprintf(fid, '%-4i%6s%8i%8s%8s%8i\t ; % 8.6f\t% 8.6f\r\n', Atoms_data{i,:});
end

fprintf(fid, '\r\n');
fprintf(fid, '[ bonds ] \r\n');
fprintf(fid, '%s\r\n','; i    j    type');

count_b = 1;
bondtype=1; % Gromacs bond type. 1 means harmonic bond, k(r-ro)^2, see manual.
% explicit_bonds = 0;
while count_b <= size(Bond_index,1)
    if explicit_bonds == 1
        if sum(ismember(Bond_index(count_b,1:2),ind_H))>0
            if strncmpi(ffname,'interface',5)
                r=0.09290;
                kb=414216;
            else
                r=0.1;
                kb=463700;
            end
        else
            if strncmpi(ffname,'interface',5)
                r=Bond_index(count_b,3)/10*1.05;
                kb=359824;
            else
                r=Bond_index(count_b,3)/10;
                kb=360000;
            end
        end
        % Normal
        Bond_order(count_b,:)= {Bond_index(count_b,1), Bond_index(count_b,2), bondtype, r, kb, ';',strtrim(char([atom(Bond_index(count_b,1)).fftype])), strtrim(char([atom(Bond_index(count_b,2)).fftype]))};
        fprintf(fid, '%-5i\t%-5i\t%-5i\t%-8.4f\t%-8.4f\t%s\t%s-%s\r\n', Bond_order{count_b,:});
        
        % Custom
        %                 Bond_order(count_b,:)= {Bond_index(count_b,1), Bond_index(count_b,2), 10, r*.95, r*1.05, r*1.05+.01 , kb, ';',strtrim(char([atom(Bond_index(count_b,1)).type])), strtrim(char([atom(Bond_index(count_b,2)).type]))};
        %                 fprintf(fid, '%-5i\t%-5i\t%-5i\t%-8.4f\t%-8.4f\t%-8.4f\t%-8.4f\t%s\t%s-%s\r\n', Bond_order{count_b,:});
        %                 fprintf(fid, '%-5i %-5i %-5i %-8.4f %-8.4f %s %s-%s\r\n', Bond_order{count_b,:});
        count_b = count_b + 1;
    else
        Bond_order(count_b,:)= {Bond_index(count_b,1), Bond_index(count_b,2), bondtype, ';', Bond_index(count_b,3)/10, strtrim(char([atom(Bond_index(count_b,1)).fftype])), strtrim(char([atom(Bond_index(count_b,2)).fftype]))};
        fprintf(fid, '%-5i %-5i %-5i %s %-8.4f %s-%s \r\n', Bond_order{count_b,:});
        count_b = count_b + 1;
    end
end

try
    if numel(Bond_order)>0
        assignin('caller','Bond_order',Bond_order);
        disp('These atom types has bonds')
        unique(Bond_order(:,end-1:end))
    end
catch
    disp('No bonds?')
end

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '[ angles ] \r\n');
fprintf(fid, '%s\r\n','; i    j   k   type');

count_a = 1;%explicit_angles = 0;
angletype=1; Angle_order={};
Angle_index=sortrows(Angle_index);
while count_a <= length(Angle_index) %nAngles;
    if explicit_angles == 1
        if sum(ismember(Angle_index(count_a,1:3),ind_H))==1
            if sum(ismember(Angle_index(count_a,1:3),ind_Mgo))>0 % Pouvreau,? Jeffery A. Greathouse,? Randall T. Cygan,? and Andrey G. Kalinichev 2017
                adeg=110;
                ktheta=50.208;
            elseif sum(ismember(Angle_index(count_a,1:3),ind_Al))>0 % Pouvreau,? Jeffery A. Greathouse,? Randall T. Cygan,? and Andrey G. Kalinichev 2017
                adeg=110; %
                ktheta=125.52; %
            else % Else orig Clayff, 2004
                adeg=110; % From most recent CHARMM prm file 116.2;
                ktheta=251.04; % since 45*4.184*2;% earlier 96.232*10; %
            end
        elseif sum(ismember(Angle_index(count_a,1:3),ind_H))==2  %             && sum(ismember(Angle_index(count_a,1:3),ind_Oh))>0 && sum(ismember(Angle_index(count_a,1:3),ind_Oct))>0
            adeg=109.47; % SPC water
            ktheta=383; % SPC water
        else % Orig Interface 2005
            adeg=Angle_index(count_a,4);
            ktheta=1422.56;
        end
        Angle_order(count_a,:)= {Angle_index(count_a,1), Angle_index(count_a,2), Angle_index(count_a,3), angletype, round(adeg,2),	ktheta, ';', strtrim(char([atom(Angle_index(count_a,1)).type])), strtrim(char([atom(Angle_index(count_a,2)).type])), strtrim(char([atom(Angle_index(count_a,3)).type]))};
        fprintf(fid, '%-5i %-5i %-5i %-5i %-6.2f %-8.4f %s %s-%s-%s\r\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    else
        Angle_order(count_a,:)= {Angle_index(count_a,1), Angle_index(count_a,2), Angle_index(count_a,3), angletype, ';', round(Angle_index(count_a,4),2), strtrim(char([atom(Angle_index(count_a,1)).fftype])), strtrim(char([atom(Angle_index(count_a,2)).fftype])), strtrim(char([atom(Angle_index(count_a,3)).fftype]))};
        fprintf(fid, '%-5i %-5i %-5i %-5i %s %-6.2f %s-%s-%s\r\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    end
end
fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

if numel(Angle_order)>0
    assignin('caller','Angle_order',Angle_order);
    disp('These atom types has angles')
    unique(Angle_order(:,end-2:end))
end

if exist('Total_charge','var')
    disp('Total charge for the .itp file was')
    round(Total_charge,5)
end

% Defining [ exclusions ]
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
% 
% if strncmpi(ffname,'clayff',5)
%     fprintf(fid, '#ifdef POSRES_Clayff \r\n');
%     fprintf(fid, '[ position_restraints ] \r\n');
%     fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
%     for i = 1:nAtoms
%         pos_res(i,:) = {num2str(i), '1', '1000', '1000', '1000'};
%         fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
%         fprintf(fid, '\r\n');
%     end
%     fprintf(fid, '#endif \r\n');
% else
%     fprintf(fid, '#ifdef POSRES_interface \r\n');
%     fprintf(fid, '[ position_restraints ] \r\n');
%     fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
%     for i = 1:nAtoms
%         pos_res(i,:) = {num2str(i), '1', '1000', '1000', '1000'};
%         fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
%         fprintf(fid, '\r\n');
%     end
%     fprintf(fid, '#endif \r\n');
% end

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

fprintf(fid, '#ifdef POSRES \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms
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
for i = 1:nAtoms
    if strncmpi([atom(i).type],'H',1)==0
        pos_res(i,:) = {num2str(i), '1', '500', '500', '500'};
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
for i = 1:nAtoms
    if ismember(i,ind_Oct)
        pos_res(i,:) = {num2str(i), '1', '500', '500', '500'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
        fprintf(fid, '\r\n');
    end
end
fprintf(fid, '#endif \r\n');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

fprintf(fid, '#ifdef POSRES_XY \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms
    if ismember(i,ind_Oct)
        if strcmp([atom(i).type],'Al') > 0 || strcmp([atom(i).type],'Mgo') > 0
            pos_res(i,:) = {num2str(i), '1', '500', '500','500'};
            fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
            fprintf(fid, '\r\n');
        end
    end
end
fprintf(fid, '#endif \r\n');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

fprintf(fid, '#ifdef POSRES_Y_10 \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms
    if ismember(i,ind_Oct)
        pos_res(i,:) = {num2str(i), '1', '1000', '10', '1000'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
        fprintf(fid, '\r\n');
    end
end
fprintf(fid, '#endif \r\n');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

fprintf(fid, '#ifdef POSRES_Y_100 \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms
    if ismember(i,ind_Oct)
        pos_res(i,:) = {num2str(i), '1', '1000', '100', '1000'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
        fprintf(fid, '\r\n');
    end
end
fprintf(fid, '#endif \r\n');

fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','[ system ]');
fprintf(fid, '%s\r\n','; name');
fprintf(fid, '%s\r\n','A mineral simulation');
fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','[ molecules ]');
fprintf(fid, '%s\r\n','; Compound        #mols');
fprintf(fid, '%s\r\n','MIN 1');
fprintf(fid, '\r\n');

fclose(fid);

[atom(strcmp([atom.type],{'Ow'})).type]=deal({'OW'});
[atom(strcmp([atom.type],{'Hw'})).type]=deal({'HW'});

atom_itp=atom;
assignin('caller','atom_itp',atom_itp);
assignin('caller','Bond_index',Bond_index);
assignin('caller','Angle_index',Angle_index);
assignin('caller','nBonds',nBonds);
assignin('caller','nAngles',nAngles);
% assignin('caller','atom',atom);
% atom = charge_atom(atom,Box_dim,ffname,watermodel,'more');
% assignin('caller','Total_charge',Total_charge);
% assignin('caller','Masses',Masses);

% toc
