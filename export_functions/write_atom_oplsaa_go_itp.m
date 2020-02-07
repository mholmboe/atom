%% write_atom_oplsaa_go_itp
% * This custom made script creates and prints a gromacs .itp file for 
% * graphene oxide using some specific OPLS/aa atom types

function write_atom_oplsaa_go_itp(atom,Box_dim,filename,varargin)

format long;
nAtoms=size(atom,2);

%
if regexp(filename,'.itp') ~= false;
    filename = filename;
else
    filename = strcat(filename,'.itp');
end

%
if nargin > 3
    short_r=cell2mat(varargin(1));
    long_r=cell2mat(varargin(2));
else
    short_r=1.25;
    long_r=1.25; % long=short since clayff
end

%
if nargin>5;
    ffname=varargin(3);
    if nargin>6;
        watermodel=varargin(4)
    else
        disp('Unknown watermodel, will try SPC/E')
        watermodel='SPC/E';
    end
    
    if strcmpi(ffname,'oplsaa_go');
        oplsaa_go_param(sort(unique([atom.type])),watermodel);
        atom = charge_opls_go_atom(atom,Box_dim,{'H' 'Oe' 'Oh'},[0.418 -0.4 -0.683],'tweak')
        Total_charge
        nrexcl=2; % See the gromacs manual
        explicit_bonds = 1;
        explicit_angles = 1;
    end
else
    disp('Unknown forcefield')
    pause
end

% [atom.temptype]=atom.type;
% atom = element_atom(atom);
% [atom.element]=atom.type;
% [atom.type]=atom.fftype
% [atom.fftype]=atom.temptype;

%
%
% atom = charge_atom(atom,Box_dim,ffname,watermodel,'more');
% assignin('caller','Total_charge',Total_charge);
atom = bond_angle_atom(atom,Box_dim,short_r,long_r);
assignin('caller','Bond_index',Bond_index);
assignin('caller','Angle_index',Angle_index);
assignin('caller','nBonds',nBonds);
assignin('caller','nAngles',nAngles);
assignin('caller','atom',atom);

%
file_title = 'Gromacs awesome itp file'; % Header in output file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

molecule_name = filename; % molecule name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Atom_label = unique([atom.type]);
%

ind_H=find(strncmp([atom.type],{'H'},1));
ind_Oh=find(strncmp([atom.type],{'Oh'},2));
ind_Oe=find(strncmp([atom.type],{'Oe'},2));
ind_Cen=find(strncmp([atom.type],{'Cen'},3));
ind_Ce=find(strncmp([atom.type],{'Ce'},2));
ind_Coh=find(strncmp([atom.type],{'Coh'},3));

fid = fopen(molecule_name, 'wt'); % open a text file that we can write into

fprintf(fid, '%s % s\r\n',';',file_title);
fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','[ moleculetype ]');
fprintf(fid, '%s % s\r\n',';','molname   nrexcl');
fprintf(fid, '%s       %d\r\n',char([atom(1).resname]),nrexcl);
fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','[ atoms ]');
fprintf(fid, '%s\r\n','; id   attype  resnr resname  atname   cgnr charge      mass');

%%
Atom_label_ID=ones(size(atom,2),1);
for i = 1:nAtoms;
    if sum(ismember(Atom_label,[atom(i).type])) > 0;
        Atom_label_ID(i,1)=find(ismember(Atom_label,[atom(i).type])==1);
    end
    %Atoms_data(i,:) = {i, char([atom(i).fftype]),1,char([atom(i).resname]),char([atom(i).type{1}(1)]),i, [atom(i).charge], Masses(Atom_label_ID(i,1))};
    Atoms_data(i,:) = {i, char([atom(i).fftype]),1,char([atom(i).resname]),char([atom(i).type]),i, [atom(i).charge], Masses(Atom_label_ID(i,1))};
    fprintf(fid, '%-4i%6s%8i%8s%8s%8i\t%8.5f\t%8.6f\r\n', Atoms_data{i,:});
end

fprintf(fid, '\r\n');
fprintf(fid, '[ bonds ] \r\n');
fprintf(fid, '%s\r\n','; i    j    type');

count_b = 1;
bondtype=1; % Gromacs bond type. 1 means harmonic bond, k(r-ro)^2, see manual.
% explicit_bonds = 0;
while count_b <= nBonds;
    if explicit_bonds == 1;
        if sum(ismember(Bond_index(count_b,1:2),ind_H))>0
            r=0.1;
            kb=462750.4;
        else
            r=Bond_index(count_b,3)/10;
            if r > 1.5
                r=1.41
            end
            kb=265265.6;
        end
        Bond_order(count_b,:)= {Bond_index(count_b,1), Bond_index(count_b,2), bondtype, r, kb, ';',strtrim(char([atom(Bond_index(count_b,1)).type])), strtrim(char([atom(Bond_index(count_b,2)).type]))};
        fprintf(fid, '%-5i %-5i %-5i %-8.4f %-8.4f %s %s-%s\r\n', Bond_order{count_b,:});
        count_b = count_b + 1;
    else
        Bond_order(count_b,:)= {Bond_index(count_b,1), Bond_index(count_b,2), bondtype, ';', Bond_index(count_b,3)/10, strtrim(char([atom(Bond_index(count_b,1)).type])), strtrim(char([atom(Bond_index(count_b,2)).type]))};
        fprintf(fid, '%-5i %-5i %-5i %s %-8.4f %s-%s \r\n', Bond_order{count_b,:});
        count_b = count_b + 1;
    end
end


if numel(Bond_order)>0;
    assignin('caller','Bond_order',Bond_order);
    disp('These atom types has bonds')
    unique(Bond_order(:,end-1:end))
end

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '[ angles ] \r\n');
fprintf(fid, '%s\r\n','; i    j   k   type');

% if strncmpi(ffname,'clayff',5);
%     disp('Removing angles with Al')
%     [Al_row,Al_col]=ind2sub(size(Angle_index),find(ismember(Angle_index,ind_Al)));
%     Angle_index(Al_row,:)=[];
%     %     %% To remove angles with 'Si'
%     %     Si_ind=find(strcmp(XYZ_labels(:,1),'Si'));
%     %     [Si_row,Si_col]=ind2sub(size(Angle_index),find(ismember(Angle_index(:,2),Si_ind)));
%     %     Angle_index(Si_row,:)=[];
% end

count_a = 1;%explicit_angles = 0;
angletype=1; Angle_order={};
Angle_index=sortrows(Angle_index);
while count_a <= length(Angle_index); %nAngles;
    if explicit_angles == 1;
        if sum(ismember(Angle_index(count_a,1:3),ind_H))>0 && sum(ismember(Angle_index(count_a,1:3),ind_Oh))>0;
            adeg=109.500; %
            ktheta=418.400; % since 45*4.184*2;% earlier 96.232*10; %
        elseif sum(ismember(Angle_index(count_a,1:3),ind_Oh))>0;
            adeg=Angle_index(count_a,4);%109.500; %
            ktheta=585.760;
        elseif sum(ismember(Angle_index(count_a,1:3),ind_Oe))>0;
            adeg=Angle_index(count_a,4);%120.00; %
            ktheta=502.080;
        else
            adeg=Angle_index(count_a,4);%120.00; %
            ktheta=527.184;
        end
        Angle_order(count_a,:)= {Angle_index(count_a,1), Angle_index(count_a,2), Angle_index(count_a,3), angletype, adeg,	ktheta, ';', strtrim(char([atom(Angle_index(count_a,1)).type])), strtrim(char([atom(Angle_index(count_a,2)).type])), strtrim(char([atom(Angle_index(count_a,3)).type]))};
        fprintf(fid, '%-5i %-5i %-5i %-5i %-8.4f %-8.4f %s %s-%s-%s\r\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    else
        Angle_order(count_a,:)= {Angle_index(count_a,1), Angle_index(count_a,2), Angle_index(count_a,3), angletype, ';', Angle_index(count_a,4), strtrim(char([atom(Angle_index(count_a,1)).type])), strtrim(char([atom(Angle_index(count_a,2)).type])), strtrim(char([atom(Angle_index(count_a,3)).type]))};
        fprintf(fid, '%-5i %-5i %-5i %-5i %s %-8.4f %s-%s-%s\r\n', Angle_order{count_a,:});
        count_a = count_a + 1;
    end
end
fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

if numel(Angle_order)>0;
    assignin('caller','Angle_order',Angle_order);
    disp('These atom types has angles')
    unique(Angle_order(:,end-2:end))
end

disp('Total charge is')
Total_charge

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

%%

if strncmpi(ffname,'oplsaa_go',5);
    fprintf(fid, '#ifdef POSRES_GO \r\n');
    fprintf(fid, '[ position_restraints ] \r\n');
    fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
    for i = 1:nAtoms;
        pos_res(i,:) = {num2str(i), '1', '1000', '1000', '1000'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
        fprintf(fid, '\r\n');
    end
    fprintf(fid, '#endif \r\n');
end

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

fprintf(fid, '#ifdef POSRES_C \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms;
    if strncmpi([atom(i).type],'C',1)==1;
    pos_res(i,:) = {num2str(i), '1', '1000', '1000', '1000'};
    fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
    fprintf(fid, '\r\n');
    end
end
fprintf(fid, '#endif \r\n');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');


fprintf(fid, '#ifdef POSRES \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms;
    pos_res(i,:) = {num2str(i), '1', '10000', '10000', '10000'};
    fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
    fprintf(fid, '\r\n');
end
fprintf(fid, '#endif \r\n');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

fprintf(fid, '#ifdef POSRES_noH \r\n');
fprintf(fid, '[ position_restraints ] \r\n');
fprintf(fid, '%s\r\n','; atom  type      fx      fy      fz');
for i = 1:nAtoms;
    if strncmpi([atom(i).type],'H',1)==0;
        pos_res(i,:) = {num2str(i), '1', '1000', '1000', '1000'};
        fprintf(fid, '%6s\t%6s\t%6s\t%6s\t%6s%\r\n', pos_res{i,:});
        fprintf(fid, '\r\n');
    end
end
fprintf(fid, '#endif \r\n');

fclose(fid);

[atom(strcmp([atom.type],{'Ow'})).type]=deal({'OW'});
[atom(strcmp([atom.type],{'Hw'})).type]=deal({'HW'});

