%% write_atom_itp.m
% * This script creates and prints a gromacs .itp file
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
% # write_atom_itp(atom,Box_dim,filename) % Basic input arguments
% # write_atom_itp(atom,Box_dim,filename,1.25,2.45,'minff','opc3') % To set cut off limits for H-bonds and all angles.

function write_minff_itp(atom,Box_dim,filename,varargin)

KANGLE=500;

format long
nAtoms=size(atom,2);

if regexp(filename,'.itp') ~= false
    filename = filename;
else
    filename = strcat(filename,'.itp');
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
round(Total_charge,5)
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
file_title = 'Gromacs awesome itp file'; % Header in output file
molecule_name = char([atom(1).resname]); % molecule name
Atom_label = unique([atom.type]);

fid = fopen(filename, 'wt'); % open a text file that we can write into

fprintf(fid, '%s % s\n',';',file_title);
fprintf(fid, '%s % s\n',';','File written by MHolmboe (michael.holmboe@umu.se)');
fprintf(fid, '\n');
fprintf(fid, '%s\n','[ moleculetype ]');
fprintf(fid, '%s % s\n',';','molname   nrexcl');
% fprintf(fid, '%s       %d\n',strrep(molecule_name(1:3),'.itp',''),nrexcl);
fprintf(fid, '%s       %d\n',molecule_name(1:3),nrexcl);
fprintf(fid, '\n');
fprintf(fid, '%s\n','[ atoms ]');
fprintf(fid, '%s\n','; id   attype  resnr resname  atname   cgnr  charge      mass');

Atom_label_ID=ones(size(atom,2),1);
for i = 1:nAtoms
    if sum(ismember(Atom_label,[atom(i).type])) > 0
        Atom_label_ID(i,1)=find(ismember(Atom_label,[atom(i).type])==1);
    end

    if isfield(atom,'mass')
        Atoms_data(i,:) = {i, char([atom(i).fftype]),[atom(i).molid],molecule_name(1:3),char([atom(i).type]),i, round([atom(i).charge],6),[atom(i).mass]};
    end
    fprintf(fid, '%-4i%6s%8i%8s%8s%8i\t% 8.6f\t% 8.6f\n', Atoms_data{i,:});
end

fprintf(fid, '\n');
fprintf(fid, '[ bonds ] \n');
fprintf(fid, '%s\n','; i    j    type');

count_b = 1;
bondtype=1; % Gromacs bond type. 1 means harmonic bond, k(r-ro)^2, see manual.
% explicit_bonds = 0;
while count_b <= size(Bond_index,1)
    if explicit_bonds == 1
        if sum(ismember(Bond_index(count_b,1:2),ind_H))>0
            r=0.09572; % 0.09789;
            kb=441050; % 463700;
        else
            r=Bond_index(count_b,3)/10;
            kb=0;
        end
        % Normal
        Bond_order(count_b,:)= {Bond_index(count_b,1), Bond_index(count_b,2), bondtype, r, kb, ';',strtrim(char([atom(Bond_index(count_b,1)).fftype])), strtrim(char([atom(Bond_index(count_b,2)).fftype]))};
        fprintf(fid, '%-5i\t%-5i\t%-5i\t%-8.4f\t%-8.4f\t%s\t%s-%s\n', Bond_order{count_b,:});

        % Custom
        %                 Bond_order(count_b,:)= {Bond_index(count_b,1), Bond_index(count_b,2), 10, r*.95, r*1.05, r*1.05+.01 , kb, ';',strtrim(char([atom(Bond_index(count_b,1)).type])), strtrim(char([atom(Bond_index(count_b,2)).type]))};
        %                 fprintf(fid, '%-5i\t%-5i\t%-5i\t%-8.4f\t%-8.4f\t%-8.4f\t%-8.4f\t%s\t%s-%s\n', Bond_order{count_b,:});
        %                 fprintf(fid, '%-5i %-5i %-5i %-8.4f %-8.4f %s %s-%s\n', Bond_order{count_b,:});
        count_b = count_b + 1;
    else
        Bond_order(count_b,:)= {Bond_index(count_b,1), Bond_index(count_b,2), bondtype, ';', Bond_index(count_b,3)/10, strtrim(char([atom(Bond_index(count_b,1)).fftype])), strtrim(char([atom(Bond_index(count_b,2)).fftype]))};
        fprintf(fid, '%-5i %-5i %-5i %s %-8.4f %s-%s \n', Bond_order{count_b,:});
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

fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '[ angles ] \n');
fprintf(fid, '%s\n','; i    j   k   type');

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
            else % Maxime Pouvreau, et al., 2019, before orig Clayff, 2004
                adeg=110; %
                ktheta=125.52; % 251.04; % since 15*4.184*2;% earlier 96.232*10; %
            end

        elseif sum(ismember(Angle_index(count_a,1:3),ind_H))==2  %             && sum(ismember(Angle_index(count_a,1:3),ind_Oh))>0 && sum(ismember(Angle_index(count_a,1:3),ind_Oct))>0
            adeg=109.47; % SPC water
            ktheta=383; % SPC water
            % elseif sum(ismember(Angle_index(count_a,1:3),ind_Al))==1 && sum(ismember(Angle_index(count_a,1:3),ind_Si))==1
            %     adeg=130; % Perfect octahedra
            %     ktheta=2000.00; % 1422.56;
            % elseif sum(ismember(Angle_index(count_a,1:3),ind_Al))>0
            %     adeg=90; % Perfect octahedra
            %     ktheta=2000.00; % 1422.56;
            % elseif sum(ismember(Angle_index(count_a,1:3),ind_Si))==1
            %     adeg=109.47; % Perfect tetrahedra
            %     ktheta=2000.00; % 1422.56;
        else % Orig Interface 2005
            adeg=Angle_index(count_a,4);
            ktheta=KANGLE;

        end
        Angle_order(count_a,:)= {Angle_index(count_a,1), Angle_index(count_a,2), Angle_index(count_a,3), angletype, round(adeg,2),	ktheta, ';', strtrim(char([atom(Angle_index(count_a,1)).type])), strtrim(char([atom(Angle_index(count_a,2)).type])), strtrim(char([atom(Angle_index(count_a,3)).type]))};
        fprintf(fid, '%-5i %-5i %-5i %-5i %-6.2f   %s %s %s-%s-%s\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    else
        Angle_order(count_a,:)= {Angle_index(count_a,1), Angle_index(count_a,2), Angle_index(count_a,3), angletype, ';', round(Angle_index(count_a,4),2), strtrim(char([atom(Angle_index(count_a,1)).fftype])), strtrim(char([atom(Angle_index(count_a,2)).fftype])), strtrim(char([atom(Angle_index(count_a,3)).fftype]))};
        fprintf(fid, '%-5i %-5i %-5i %-5i %s %-6.2f %s-%s-%s\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    end
end
fprintf(fid, '\n');
fprintf(fid, '\n');

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
%         fprintf(fid, '\n');
%         fprintf(fid, '\n');
%         fprintf(fid, '[ exclusions ] \n');
%         fprintf(fid, '%s\n','; i    j   k   type');
%
%         count_excl = 1;
%         Excl_index=[Angle_index(:,1) Angle_index(:,2) Angle_index(:,3); Angle_index(:,2) Angle_index(:,3) Angle_index(:,1); Angle_index(:,2) Angle_index(:,3) Angle_index(:,1)];
%         while count_excl <= length(Excl_index);
%             Excl_order(count_excl,:)= {Excl_index(count_excl,1), Excl_index(count_excl,2), Excl_index(count_excl,3),';', strtrim(char(XYZ_labels(Excl_index(count_excl,1)))), strtrim(char(XYZ_labels(Excl_index(count_excl,2)))), strtrim(char(XYZ_labels(Excl_index(count_excl,3))))};
%             fprintf(fid, '%-5i %-5i %-5i %s %s-%s-%s\n', Excl_order{count_excl,:});
%             count_excl = count_excl + 1;
%         end
%
%     end

%%%%%%%%%%%%%%%%%%

fprintf(fid, '#ifdef POSRES_XY_MMT_1  \n');
fprintf(fid, '[ position_restraints ] \n');
fprintf(fid, '%s\n','; atom  type      fx      fy      fz');

for i = 1:nAtoms
    if ismember(i,ind_Oct)
        pos_res(i,:) = {num2str(i), '1', '100', '100', '10000'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\n', pos_res{i,:});
        fprintf(fid, '\n');
    end
end
fprintf(fid, '#endif \n');

fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, '#ifdef POSRES_XY_MMT_2  \n');
fprintf(fid, '[ position_restraints ] \n');
fprintf(fid, '%s\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms
    if ismember(i,ind_Oct)
        pos_res(i,:) = {num2str(i), '1', '1000', '1000', '0'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\n', pos_res{i,:});
        fprintf(fid, '\n');
    end
end
fprintf(fid, '#endif \n');

fprintf(fid, '\n');
fprintf(fid, '\n');

%%%%%%%%%%%%

if strncmpi(ffname,'minff',5)
    fprintf(fid, '#ifdef POSRES_MINFF \n');
    fprintf(fid, '[ position_restraints ] \n');
    fprintf(fid, '%s\n','; atom  type      fx      fy      fz');
    for i = 1:nAtoms
        pos_res(i,:) = {num2str(i), '1', '1000', '1000', '1000'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\n', pos_res{i,:});
        fprintf(fid, '\n');
    end
    fprintf(fid, '#endif \n');
end

fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, '#ifdef POSRES \n');
fprintf(fid, '[ position_restraints ] \n');
fprintf(fid, '%s\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms
    pos_res(i,:) = {num2str(i), '1', '1000', '1000', '1000'};
    fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\n', pos_res{i,:});
    fprintf(fid, '\n');
end
fprintf(fid, '#endif \n');

fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, '#ifdef POSRES_noH \n');
fprintf(fid, '[ position_restraints ] \n');
fprintf(fid, '%s\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms
    if strncmpi([atom(i).type],'H',1)==0
        pos_res(i,:) = {num2str(i), '1', '1000', '1000', '1000'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\n', pos_res{i,:});
        fprintf(fid, '\n');
    end
end
fprintf(fid, '#endif \n');

fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, '#ifdef POSRES_XYZ \n');
fprintf(fid, '[ position_restraints ] \n');
fprintf(fid, '%s\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms
    if strncmpi([atom(i).type],'H',1)==0
        pos_res(i,:) = {num2str(i), '1', '500', '500', '500'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\n', pos_res{i,:});
        fprintf(fid, '\n');
    end
end
fprintf(fid, '#endif \n');

fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, '#ifdef POSRES_XY \n');
fprintf(fid, '[ position_restraints ] \n');
fprintf(fid, '%s\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms
    if strncmpi([atom(i).type],'H',1)==0
        %         if strcmp([atom(i).type],'Al') > 0 || strcmp([atom(i).type],'Mgo') > 0
        pos_res(i,:) = {num2str(i), '1', '1000', '1000','0'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\n', pos_res{i,:});
        fprintf(fid, '\n');
        %         end
    end
end
fprintf(fid, '#endif \n');

fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, '#ifdef POSRES_Oct_500 \n');
fprintf(fid, '[ position_restraints ] \n');
fprintf(fid, '%s\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms
    if ismember(i,ind_Al)
        pos_res(i,:) = {num2str(i), '1', '500', '500', '500'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\n', pos_res{i,:});
        fprintf(fid, '\n');
    end
end
fprintf(fid, '#endif \n');

fprintf(fid, '\n');
fprintf(fid, '\n');

fprintf(fid, '#ifdef POSRES_Y_100 \n');
fprintf(fid, '[ position_restraints ] \n');
fprintf(fid, '%s\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms
    if ismember(i,ind_Oct)
        pos_res(i,:) = {num2str(i), '1', '1000', '100', '1000'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\n', pos_res{i,:});
        fprintf(fid, '\n');
    end
end
fprintf(fid, '#endif \n');

fprintf(fid, '\n');
fprintf(fid, '\n');


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