%% add2atom.m
% This function appends so-called XYZ atomtype labels and XYZ data to an existing 
% atom struct
% 
% * XYZ_labels is a nAtoms x 1 (cell) vector of atomtype names
% * XYZ_data is a nAtoms x 3 data matrix holding the x,y,z coordinates
% * resname is the optional Resname you choose
% * in_atom is the optional (requires resname) existing atom struct for which 
% the XYZ data should be appended to
% 
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = add2atom(XYZ_labels,XYZ_data) % Basic input arguments
% # atom = add2atom(XYZ_labels,XYZ_data,'MOL') % Will set the residue name to MOL
% # atom = add2atom(XYZ_labels,XYZ_data,'LAC',in_atom) % Will set the residue name to LAC and the MolID to whatever comes after the MolID in in_atom
%
function atom = add2atom(XYZ_labels,XYZ_data,varargin)

if numel(XYZ_data)==1
    if iscell(XYZ_labels)==0;XYZ_labels={XYZ_labels}; end
end
if numel(XYZ_data)==1
    XYZ_data(1)=XYZ_data(1);
    XYZ_data(2)=XYZ_data(1);
    XYZ_data(3)=XYZ_data(1);
end
if nargin > 2
    resname=varargin{1};
    if nargin > 3
        in_atom=varargin{2};
    else
        in_atom=[];
    end
else
    resname=XYZ_labels;
    in_atom=[];
end
nAtoms=length(XYZ_data);
nmol=1; MolID=zeros(nAtoms,1);
if strncmpi(XYZ_labels(1),'Ow',2)
    MolID(1:3:end)=ceil((nmol/3):3:nAtoms)';
    MolID(2:3:end)=ceil((nmol/3):3:nAtoms)';
    MolID(3:3:end)=ceil((nmol/3):3:nAtoms)';
else
    MolID=nmol:1:nAtoms;
end
if isfield(in_atom,'molid')
    MolID=MolID+in_atom(end).molid;
    index=in_atom(end).index;
else
    in_atom=[];
    index=0;
end
if isfield(in_atom,'molid') == false && sum(size(unique(XYZ_labels),1)) > 3
    MolID(:)=1;
    index=0;
elseif isfield(in_atom,'molid') == true && sum(size(unique(XYZ_labels),1)) > 3
    MolID(:)=1+in_atom(end).molid;
    index=in_atom(end).index;
end
nmol=MolID(1);
% Put everything in the atom struct
first_in=[1];last_in=[];
for i=1:size(XYZ_data,1)
    if i > 1 && MolID(i) ~= MolID(i-1)
        nmol=nmol+1;
        atom(i).molid=nmol;
        first_in(atom(i).molid,1)=i; last_in(atom(i).molid-1,1)=i-1;
    elseif i > 1
        atom(i).molid=atom(i-1).molid;
    elseif i == 1
        atom(i).molid=MolID(1);
    end
    if strcmpi(resname,'same')
        disp('Same resname as type')
        atom(i).resname=XYZ_labels(i);
    else
        atom(i).resname=resname;
    end
    atom(i).type=XYZ_labels(i,:);
    atom(i).fftype=XYZ_labels(i,:);
    atom(i).index=index+mod(i,100000);
    atom(i).neigh.type  = {};
    atom(i).neigh.index  = zeros(6,1);
    atom(i).neigh.dist  = zeros(6,1);
    atom(i).bond.type  = zeros(6,1);
    atom(i).bond.index  = zeros(6,1);
    atom(i).angle.type  = zeros(6,1);
    atom(i).angle.index  = zeros(6,1);
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
atom=[in_atom atom];
assignin('caller','XYZ_data',[[atom.x]' [atom.y]' [atom.z]']);
assignin('caller','XYZ_labels',[atom.type]');
assignin('caller','nAtoms',nAtoms)
assignin('caller','MolID',MolID)
disp('add2atom done!')