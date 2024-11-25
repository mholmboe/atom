%% hbond_atom.m
% * This function tries to calculate the number of hydrogen bonds
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=hbond_atom(atom,Box_dim) % Basic input arguments
% # atom=hbond_atom(atom,Box_dim,hb_dist) % Sets the X-H-X distance
% # atom=hbond_atom(atom,Box_dim,hb_dist,hbangle) % Sets the X-H-X angle

function atom=hbond_atom(atom,Box_dim,varargin)

format short
outfilename='hbond_conf.pdb';
hb_dist=3.5; % Å
hb_angle=30; % max allowed deg between H-D---A
DAgroup='O'; % Only one Donor/acceptor group at the moment, otherwise add indexes to ind1 (which must contain indexes of both the D's and A's)
DHgroup='H'; % The donor-H's atomtypes, like 'HW'

ind1=find(strncmpi([atom.type],DHgroup,1)); % indexes of the acceptors
ind2=find(strncmpi([atom.type],DAgroup,1)); % indexes of the donor-H's
ind2=unique(sort([ind2 find(strncmpi([atom.type],'Oh',2))])); % indexes of the donor-H's]

dist_step=0.1;
angle_step=1;
% bins_HBD=zeros(1,size(atom,2));
bins_DIST=zeros(1,length(0:dist_step:hb_dist)-1);
bins_ANGLE=zeros(1,length(0:angle_step:hb_angle)-1);

atom1=atom(ind1); % Extract the solute atomtype
atom2=atom(ind2); % Extract the ligands atomtypes

filter=[];
if nargin>2
    filter=varargin{1};
end

if numel(filter)>0
    % 'Filter' the structure here, with respect to donor-H's...
    % Select Hw close to ind1 atoms
    dist_matrix = dist_matrixes_atom(atom1,atom2,Box_dim);
    sel_ind_matrix=dist_matrix < hb_dist;
    sel_ind=find(sum(sel_ind_matrix)>0);
    sel_ind=intersect(1:3:size(atom2,2),sel_ind);
    sel_ind=sort([sel_ind sel_ind+1 sel_ind+2]);
    % % Invert the selection
    % sel_ind=setdiff(1:size(atom2,2),sel_ind);
    atom2=atom2(sel_ind);
end

Zslice=[];
if nargin>3
    Zslice=varargin{2};
end

if numel(Zslice)>1
    %% Filter Hw with certain z coordinate
    atom2=atom2([atom.z]>min(Zslice)&[atom.z]<=max(Zslice));
end


%% This is the slow step...
temp_atom = bond_angle_type(atom1,atom2,Box_dim,0.8,hb_dist-1,120,0,'min_angle');

%% Put together the Angle_index matrix
Angle_index=cell2mat(arrayfun(@(x) x.neigh.angle, temp_atom(1:end),'UniformOutput',0)');
temp_Angle_index=Angle_index;
ind_reverse=find(Angle_index(:,5)>1.2);
temp_Angle_index(ind_reverse,[1 3 5:12])=Angle_index(ind_reverse,[3 1 6 5 10 11 12 7 8 9]);
Angle_index=temp_Angle_index;

%% Check the D-H--A angles
HDA_angle=cell2mat(arrayfun(@(x) rad2deg(atan2(norm(cross([Angle_index(x,10:12)]-[Angle_index(x,7:9)],[0 0 0]-[Angle_index(x,7:9)])),...
    dot([Angle_index(x,10:12)]-[Angle_index(x,7:9)],[0 0 0]-[Angle_index(x,7:9)]))), 1:size(Angle_index,1),'UniformOutput',0)');
HDA_angle_ind=find(HDA_angle<hb_angle);

%% Check the D--A distances
DA_dist=cell2mat(arrayfun(@(x) pdist2(Angle_index(x,7:9),Angle_index(x,10:12)), 1:size(Angle_index,1),'UniformOutput',0)');
DA_dist_ind=find(DA_dist<hb_dist);
HDA_ind=intersect(HDA_angle_ind,DA_dist_ind);

%% Cut out the cutout's
Angle_index=Angle_index(HDA_ind,:);

%% Collect the donor Hbonds
bins_HBD(1,:)=histcounts(sort([Angle_index(:,1)]),size(atom,2))';
bins_DIST=histcounts(Angle_index(:,5),0:dist_step:hb_dist)';
bins_ANGLE=histcounts(180-Angle_index(:,4),0:angle_step:hb_angle)';

HBD.HB1=find(bins_HBD(1,:)==1);
HBD.HB2=find(bins_HBD(1,:)==2);
HBD.HB3=find(bins_HBD(1,:)==3);
HBD.HB4=find(bins_HBD(1,:)==4);
HBD.HB5=find(bins_HBD(1,:)==5);
nHD_bonds(1)=(1*length([HBD.HB1])+2*length([HBD.HB2])+3*length([HBD.HB3])...
    +4*length([HBD.HB4])+5*length([HBD.HB5]))/size(atom2,2);
disp('Total number of Donator H-bonds');
nHD_bonds(end)

%% Collect the acceptor Hbonds
bins_HBA(1,:)=histcounts(sort([Angle_index(:,3)]),size(atom,2))';
bins_DIST=histcounts(Angle_index(:,6),0:dist_step:hb_dist)';
bins_ANGLE=histcounts(180-Angle_index(:,4),0:angle_step:hb_angle)';

HBA.HB1=find(bins_HBA(1,:)==1);
HBA.HB2=find(bins_HBA(1,:)==2);
HBA.HB3=find(bins_HBA(1,:)==3);
HBA.HB4=find(bins_HBA(1,:)==4);
HBA.HB5=find(bins_HBA(1,:)==5);
nHA_bonds(1)=(1*length([HBA.HB1])+2*length([HBA.HB2])+3*length([HBA.HB3])...
    +4*length([HBA.HB4])+5*length([HBA.HB5]))/size(atom2,2);
disp('Total number of Acceptor H-bonds');
nHA_bonds(end)

disp('Total number of H-bonds');
nHD_bonds(end)+nHA_bonds(end)


% % % Another way of counting the number of HB's
% % nH_bondsv2=size(Angle_index,1)*2/size(atom1,2)

assignin('caller','bins_HB',bins_HBD);
assignin('caller','bins_DIST',bins_DIST);
assignin('caller','bins_ANGLE',bins_ANGLE);
assignin('caller','HB',HBD);
% assignin('caller','Bond_index',Bond_index);
assignin('caller','Angle_index',Angle_index);

% Print PDB and CONECT files, showing the HB's
%% Write a pdb without the CONECT section
write_atom_pdb(atom,Box_dim,outfilename)

%% Print the CONECT section with th HB bonds separately
HB_Bond_index=[Angle_index(find(Angle_index(:,5)>1.2),1:2) Angle_index(find(Angle_index(:,5)>1.2),end); Angle_index(find(Angle_index(:,6)>1.2),2:3) Angle_index(find(Angle_index(:,6)>1.2),end)];
HB_Bond_index=sortrows(HB_Bond_index,1);

nAtoms=size(atom,2);
Bond_index=HB_Bond_index(:,1:2);
B=[Bond_index(:,1:2); Bond_index(:,2) Bond_index(:,1)];
b1=sortrows(B);

fid = fopen(outfilename, 'a+');
for i=1:max(b1(:,1))
    ind=find(b1(:,1)==i);
    b2=b1(ind,2);
    fprintf(fid,'CONECT%5i%5i%5i%5i%5i%5i%5i',[i;b2]);
    fprintf(fid,'\r\n');
end

fprintf(fid,'MASTER    %5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i\r\n',[0    0    0    0    0    0    0    0 nAtoms    0 i    0]);
fprintf(fid,'END');
fclose(fid);


