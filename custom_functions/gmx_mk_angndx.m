%% gmx_mk_angndx.m
% * This function helps you print custom gromacs angle .ndx files with
% triplets indexes. 
% * atomtype1-3 can either be atom names like 'Al', or an array of indexes 
% * like [1:10:960] or [1 2 3 4 5 8].
% * rmin and rmax can be used to play around with bond lengths.
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% * gmx_mk_angndx(atom,Box_dim,'Mgo','Ohmg','H','MgOH')
% * gmx_mk_angndx(atom,Box_dim,'Mgo','Ohmg','H','MgOH',1,1.95)
% * gmx_mk_angndx(atom,Box_dim,IndexArray_atom1,'Ohmg',IndexArray_atom3,'angle_index_test')

function ind = gmx_mk_angndx(atom,Box_dim,atomtype1,atomtype2,atomtype3,groupname,varargin)

if nargin > 6
    rmin=varargin{1};
    rmax=varargin{2};
else
    rmin=1.25;
    rmax=2.25;
end

orig_atom=atom;
if ischar(atomtype1) && ischar(atomtype2) && ischar(atomtype3)
    atom=atom(ismember([atom.type],{atomtype1 atomtype2 atomtype3}));
end

if ischar(atomtype1)
    ind_atom1=find(ismember([atom.type],atomtype1));
else
    ind_atom1=atomtype1;
end
if ischar(atomtype2)
    ind_centertype=find(ismember([atom.type],atomtype2));
else
    ind_centertype=atomtype2;
end
if ischar(atomtype3)
    ind_atom3=find(ismember([atom.type],atomtype3));
else
    ind_atom3=atomtype3;
end

[atom(ind_centertype).type]=deal({'Center'});
atom=bond_angle_atom(atom,Box_dim,1.25,2.25,'more');
Center_angles=Center_angles(ismember(Center_angles(:,2),ind_centertype),:);

i=1;
while i<length(Center_angles)+1
    if sum(ismember(Center_angles(i,1:3),[ind_atom1 ind_atom3]))<2
        Center_angles(i,:)=[];
    else
        i=i+1;
    end
end

angndx=Center_angles(:,1:3);
angndx=reshape(angndx',1,[]);

ind=[orig_atom([atom(angndx).index]).index];

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