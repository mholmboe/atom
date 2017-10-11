%% atom2make_ndx.m
% * This function can help you print custom gromacs .ndx files
% * Likely best used interactively... se below
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom2make_ndx(filename,groupname,atomtypes,molid)


function atom2make_ndx(filename,groupname,atomtypes,molid)

% clear all; format compact;
% You could set the group name here manually if you want
% filename='preem.gro'
% groupname='FIX_MMT_1';
% atomtypes={'Al' 'Mgo'};
% molid=1;

atom=import_atom(filename);
ind_atomtypes=ismember([atom.type],atomtypes);
ind_molid=ismember([atom.molid],molid);
ind=ind_molid+ind_atomtypes;
ind(ind<2)=0;ind(ind>0)=1;
ind=find(ind==1);

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
