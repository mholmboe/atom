%% xyz2atom.m
% * This function can be used to add XYZ data (like from a .xyz structure 
% file) to the atom struct format.
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = xyz2atom(XYZ_labels,XYZ_data,Box_dim,resname,in_atom)
%
function atom = xyz2atom(XYZ_labels,XYZ_data,Box_dim,resname,in_atom)

nAtoms=length(XYZ_data);

natom=1; molid=zeros(nAtoms,1);
if strncmpi(XYZ_labels(1),'Ow',2)
    molid(1:3:end)=ceil((natom/3):3:nAtoms)';
    molid(2:3:end)=ceil((natom/3):3:nAtoms)';
    molid(3:3:end)=ceil((natom/3):3:nAtoms)';
else
    molid=natom:1:nAtoms;
end

if isfield(in_atom,'molid')
    molid=molid+in_atom(end).molid;
    index=in_atom(end).index;
else
    in_atom=[];
    index=0;
end

if isfield(in_atom,'molid') == false && sum(size(unique(XYZ_labels),1)) > 3
    molid(:)=1;
    index=0;
elseif isfield(in_atom,'molid') == true && sum(size(unique(XYZ_labels),1)) > 3
    molid(:)=1+in_atom(end).molid;
    index=in_atom(end).index;
end

natom=molid(1);


% This seems also to work for the coordinates but is 10% slower...
% X=num2cell(new_XYZ_data(:,1));
% Y=num2cell(new_XYZ_data(:,2));
% Z=num2cell(new_XYZ_data(:,3));
% [new_atom.x]=X{:};
% [new_atom.y]=Y{:};
% [new_atom.z]=Z{:};

first_in=[1];last_in=[];
for i=1:nAtoms
    if i > 1 && molid(i) ~= molid(i-1)
        natom=natom+1;
        atom(i).molid=natom;
        first_in(atom(i).molid,1)=i; last_in(atom(i).molid-1,1)=i-1;
    elseif i > 1
        atom(i).molid=atom(i-1).molid;
    elseif i == 1
        atom(i).molid=molid(1);
    end
    atom(i).resname=resname;%XYZ_labels(i);
    atom(i).type=XYZ_labels(i);
    atom(i).fftype=XYZ_labels(i);
    atom(i).index=index+mod(i,100000);
    atom(i).neigh.type  = {};
    atom(i).neigh.index  = zeros(6,1);
    atom(i).neigh.dist  = zeros(6,1);
    atom(i).bond.type  = zeros(6,1);
    atom(i).bond.index  = zeros(6,1);
    atom(i).x=XYZ_data(i,1);
    atom(i).y=XYZ_data(i,2);
    atom(i).z=XYZ_data(i,3);
    atom(i).vx=nan;
    atom(i).vy=nan;
    atom(i).vz=nan;
end

if sum(size(unique(XYZ_labels),1)) > 3
    [atom.resname]=deal({resname});
end

atom=[in_atom, atom];

% XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
% XYZ_labels=[atom.type]';

%assignin('caller','atom',atom);
assignin('caller','XYZ_data',[[atom.x]' [atom.y]' [atom.z]']);
assignin('caller','XYZ_labels',[atom.type]');
assignin('caller','nAtoms',nAtoms)
assignin('caller','Box_dim',Box_dim)
assignin('caller','molid',molid)

disp('add2atom done!')
