%% slice_triclinic_atom.m
% * This function slices the atoms into the triclinic box
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = slice_triclinic_atom(atom,Box_dim) % Basic input arguments
% # atom = slice_triclinic_atom(atom,Box_dim,'xy') % Will only slice along the xy dimensions
%
function atom = slice_triclinic_atom(atom,Box_dim,varargin)
% Which one is fastest? Ortogonal or triclinic?

if nargin==2
    disp('will try to slice the triclinic Box_dim!!!')
    
    if size(Box_dim(1,:),2)>3
        xy=Box_dim(6); xz=Box_dim(8); yz=Box_dim(9);
    else
        xy=0;xz=0;yz=0;
    end
    orto=orto_atom(atom,Box_dim);
    
    indxlo=find([orto.xfrac]<0);
    indxhi=find([orto.xfrac]>1);
    
    indylo=find([orto.yfrac]<0);
    indyhi=find([orto.yfrac]>1);
    
    indzlo=find([orto.zfrac]<0);
    indzhi=find([orto.zfrac]>1);
    
    ind = unique([indxlo indxhi indylo indyhi indzlo indzhi]);
    orto(ind)=[];
    orto=update_atom(orto);
    atom = triclinic_atom(orto,orto_Box_dim,[xy xz yz],'tilt');
    
elseif nargin == 3
    
    disp('will try to slice the triclinic Box_dim along xy!!!')
    xy=Box_dim(6); xz=Box_dim(8); yz=Box_dim(9);
    orto=orto_atom(atom,Box_dim);
    
    indxlo=find([orto.xfrac]<0);
    indxhi=find([orto.xfrac]>1);
    
    indylo=find([orto.yfrac]<0);
    indyhi=find([orto.yfrac]>1);
    
    indzlo=find([orto.zfrac]<0);
    indzhi=find([orto.zfrac]>1);
    
    ind = unique([indxlo indxhi indylo indyhi]); % indzlo indzhi]);
    orto(ind)=[];
    orto=update_atom(orto);
    atom = triclinic_atom(orto,orto_Box_dim,[xy xz yz],'tilt');
    
end

try
    if isfield(atom,'xfrac')
        atom=rmfield(atom,'xfrac');
        atom=rmfield(atom,'yfrac');
        atom=rmfield(atom,'zfrac');
    end
catch
end

%assignin('caller','atom',atom);