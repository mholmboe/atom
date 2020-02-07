%% atom2make_ndx.m
% * This special function can help you print custom gromacs .ndx files, you
% * just have to find a way to find the indexes you want.
% * Likely best used interactively... se below
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom2make_ndx(filename,groupname,atomtypes,molid)

function atom2make_ndx(filename,groupname,atomtypes,molid)

% clear all; format compact;
% You could set the group name here manually if you want
% filename='evap_0.gro'
% groupname='MMT';
% atomtypes={'Mgo' 'Ohmg'};
% molid=[1 2 3];

atom=import_atom(filename);
ind_atomtypes=ismember([atom.type],atomtypes);
ind_molid=ismember([atom.molid],molid);
ind=ind_molid+ind_atomtypes;
ind(ind<2)=0;ind(ind>0)=1;
ind=find(ind==1);

% a 7-9 | a 27-29
% 
% name 13 low_Ob
% 
% a 17-19 | a 37-39
% 
% name 14 hi_Ob

filename='npt.ndx';
groupname='dihedrals';
ind=sort([[6:40:960] [10:40:960] [16:40:960] [20:40:960]]);
ind=sort([[1:10:5000] [2:10:5000] [3:10:5000] [4:10:5000]]);

% groupname='MOH_angle';
% indMg=find(strncmpi([atom.type],'Al',2));
% indOhmg=find(ismember([atom.type],'Oh'));
% ind=sort([indMg indOhmg indOhmg+4]);


% Format the ind vector to have 15 entries per row
ext_ind=zeros(1,15*ceil(length(ind)/15));
ext_ind(1:length(ind))=ind;
ext_ind=reshape(ext_ind,15,[])';

% Print the index file
fid = fopen(strcat(groupname,'.ndx'), 'wt');
fprintf(fid, '%s %s %s\r\n','[',groupname,']');
for i = 1:size(ext_ind,1)
    row=ext_ind(i,:);
    row(row==0)=[];
    if max(row)<1000
    fprintf(fid, '%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i\r\n', row);
    elseif max(row)<10000
    fprintf(fid, '%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i\r\n', row);
    elseif max(row)<100000
    fprintf(fid, '%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i\r\n', row);
    end
end
fprintf(fid, '\r\n');
fclose(fid);
