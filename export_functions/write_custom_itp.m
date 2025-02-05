%% write_custom_itp.m
% * This script creates and prints a custom gromacs .itp file
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_custom_itp(atom,Box_dim,filename) % Basic input arguments
% # write_custom_itp(atom,Box_dim,filename,1.25,2.25) % To set the H-O and M-O cutoff radii, resp.

function write_custom_itp(atom,Box_dim,filename,varargin)

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


disp('Forcefield not stated, will make some assumptions then...')
ffname='custom';
watermodel='SPC'; % SPC/E, depreceated
Total_charge=sum([atom.charge])
round2dec(Total_charge,5)
%         pause
nrexcl=1; % See the gromacs manual
explicit_bonds = 0;
explicit_angles = 0;
MolId=atom(1).molid;

if isfield(atom,'element')==0
    element=element_atom(atom);
    [atom.element]=element.type;
end
[atom.fftype]=atom.type;
atom = bond_angle_dihedral_atom(atom,Box_dim,maxrshort,maxrlong);
atom = mass_atom(atom);
atom=update_atom(atom);
%
file_title = 'Gromacs awesome itp file'; % Header in output file
molecule_name = char([atom(1).resname]); % molecule name
Atom_label = unique([atom.type]);

fid = fopen(filename, 'wt'); % open a text file that we can write into

% fprintf(fid, '%s % s\n',';',file_title);
fprintf(fid, '%s % s\n',';','File written by MHolmboe (michael.holmboe@umu.se)');
fprintf(fid, '\n');
fprintf(fid, '%s\n','[ moleculetype ]');
fprintf(fid, '%s % s\n',';','molname   nrexcl');
% fprintf(fid, '%s       %d\n',strrep(molecule_name(1:3),'.itp',''),nrexcl);
fprintf(fid, '%s       %d\n',molecule_name(1:3),nrexcl);
fprintf(fid, '\n');
fprintf(fid, '%s\n','[ atoms ]');
fprintf(fid, '%s\n','; id   attype  resnr resname  atname   cgnr  charge      mass');

Atom_label_ID=ones(size(atom,2),1);sum_charge=0;
for i = 1:nAtoms
    if sum(ismember(Atom_label,[atom(i).type])) > 0
        Atom_label_ID(i,1)=find(ismember(Atom_label,[atom(i).type])==1);
    end
    charge=round2dec([atom(i).charge],5)%+0.00328402;
    sum_charge=sum_charge+charge;
    Atoms_data(i,:) = {i, char([atom(i).fftype]),[atom(i).molid],molecule_name(1:3),char([atom(i).element]),i, charge,[atom(i).mass],';',sum_charge};
    % fprintf(fid, '%-4i\t%6s\t%8i\t%8s\t%8s\t%8i\t%8.6f\t%-8.6f\n', Atoms_data{i,:});
    fprintf(fid, '%6i%11s%9i%5s%7s%7i\t%8.5f\t%8.5f\t%5s\t%8.5f\n', Atoms_data{i,:});
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
            kb=360000;
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
        fprintf(fid, '%-5i %-5i %-5i %s %-8.4f %s-%s\n', Bond_order{count_b,:});
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

%% To include a generic 1-4 pairlist
if length(Pairlist)>0

    fprintf(fid, '\n');
    fprintf(fid, '[ pairs ] \n');
    fprintf(fid, '%s\n','; ai	aj	funct	c6	c12 or');
    % fprintf(fid, '%s\n','; ai	aj	funct	fudgeQQ	q1	q2	c6	c12');

    count_p = 1;%explicit_angles = 0;
    Pairlisttype=1; Pair_order={};
    while count_p <= length(Pairlist) %nAngles;
        Pair_order(count_p,:)= {Pairlist(count_p,1), Pairlist(count_p,2), Pairlisttype, ';',strtrim(char([atom(Pairlist(count_p,1)).fftype])), strtrim(char([atom(Pairlist(count_p,2)).fftype]))};
        fprintf(fid, '%-5i %-5i %-5i %s %s-%s\n', Pair_order{count_p,:});
        count_p=count_p+1;
    end
end

fprintf(fid, '\n');
fprintf(fid, '[ angles ] \n');
fprintf(fid, '%s\n','; i    j   k   type');

count_a = 1;%explicit_angles = 0;
angletype=5; Angle_order={};
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
            if strncmpi(ffname,'interface',5)
                adeg=126.00;
                ktheta=376.56;
            end

        elseif sum(ismember(Angle_index(count_a,1:3),ind_H))==2  %             && sum(ismember(Angle_index(count_a,1:3),ind_Oh))>0 && sum(ismember(Angle_index(count_a,1:3),ind_Oct))>0
            adeg=109.47; % SPC water
            ktheta=383; % SPC water
        end
        Angle_order(count_a,:)= {Angle_index(count_a,1), Angle_index(count_a,2), Angle_index(count_a,3), angletype, round2dec(adeg,2),	ktheta, ';', strtrim(char([atom(Angle_index(count_a,1)).type])), strtrim(char([atom(Angle_index(count_a,2)).type])), strtrim(char([atom(Angle_index(count_a,3)).type]))};
        fprintf(fid, '%-5i %-5i %-5i %-5i %-6.2f   %-8.2f %s %s-%s-%s\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    else
        Angle_order(count_a,:)= {Angle_index(count_a,1), Angle_index(count_a,2), Angle_index(count_a,3), angletype, ';', round2dec(Angle_index(count_a,4),2), strtrim(char([atom(Angle_index(count_a,1)).fftype])), strtrim(char([atom(Angle_index(count_a,2)).fftype])), strtrim(char([atom(Angle_index(count_a,3)).fftype]))};
        fprintf(fid, '%-5i %-5i %-5i %-5i %s %-6.2f %s-%s-%s\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    end
end
fprintf(fid, '\n');

if numel(Angle_order)>0
    assignin('caller','Angle_order',Angle_order);
    disp('These atom types has angles')
    unique(Angle_order(:,end-2:end))
end

fprintf(fid, '[ dihedrals ] \n');
fprintf(fid, '%s\n','; i    j   k   type');

count_d = 1;
dihedraltype=9; Dihedral_order={};
Dihedral_index=sortrows(Dihedral_index);
while count_d <= length(Dihedral_index)
    Dihedral_order(count_d,:)= {Dihedral_index(count_d,1), Dihedral_index(count_d,2), Dihedral_index(count_d,3), Dihedral_index(count_d,4), dihedraltype, ';',...
        strtrim(char([atom(Dihedral_index(count_d,1)).type])), strtrim(char([atom(Dihedral_index(count_d,2)).type])), strtrim(char([atom(Dihedral_index(count_d,3)).type])), strtrim(char([atom(Dihedral_index(count_d,4)).type]))};
    fprintf(fid, '%-5i %-5i %-5i %-5i %-5i %s %s-%s-%s-%s\n', Dihedral_order{count_d,:});
    count_d = count_d + 1;
end
fprintf(fid, '\n');

if numel(Dihedral_order)>0
    assignin('caller','Dihedral_order',Dihedral_order);
    disp('These atom types has dihedrals')
    unique(Dihedral_order(:,end-2:end))
end

fprintf(fid, '[ dihedrals ] \n');
fprintf(fid, '%s\n','; i    j   k   type');

count_d = 1;
dihedraltype=2; Improper_dihedral_order={};
Improper_dihedral_index=sortrows(Improper_dihedral_index);
while count_d <= length(Improper_dihedral_index)
    Improper_dihedral_order(count_d,:)= {Improper_dihedral_index(count_d,1), Improper_dihedral_index(count_d,2), Improper_dihedral_index(count_d,3), Improper_dihedral_index(count_d,4), dihedraltype, ';',...
        strtrim(char([atom(Improper_dihedral_index(count_d,1)).type])), strtrim(char([atom(Improper_dihedral_index(count_d,2)).type])), strtrim(char([atom(Improper_dihedral_index(count_d,3)).type])), strtrim(char([atom(Improper_dihedral_index(count_d,4)).type]))};
    fprintf(fid, '%-5i %-5i %-5i %-5i %-5i %s %s-%s-%s-%s\n', Improper_dihedral_order{count_d,:});
    count_d = count_d + 1;
end
fprintf(fid, '\n');

if numel(Improper_dihedral_order)>0
    assignin('caller','Improper_dihedral_order',Improper_dihedral_order);
    disp('These atom types has improper dihedrals')
    unique(Improper_dihedral_order(:,end-2:end))
end



if exist('Total_charge','var')
    disp('Total charge for the .itp file was')
    round2dec(sum(cell2mat(Atoms_data(:,7))),5)
    round2dec(Total_charge,5)
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

fclose(fid);

[atom(strcmp([atom.type],{'Ow'})).type]=deal({'OW'});
[atom(strcmp([atom.type],{'Hw'})).type]=deal({'HW'});

atom_itp=atom;
assignin('caller','Atoms_data',Atoms_data);
assignin('caller','atom_itp',atom_itp);
assignin('caller','Bond_index',Bond_index);
assignin('caller','Angle_index',Angle_index);
assignin('caller','Dihedral_index',Dihedral_index);
assignin('caller','Improper_dihedral_index',Improper_dihedral_index);
assignin('caller','nBonds',nBonds);
assignin('caller','nAngles',nAngles);
assignin('caller','nDihedrals',nDihedrals);
assignin('caller','nImpropers',nImpropers);
