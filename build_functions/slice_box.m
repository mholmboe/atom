%% slice_box.m
% * This function checks if the coordinates is within the specified limits,
% and if not sets the x,y,z to nan,nan,nan. See also slice_atom/slice_molid
% functions, which slices taking the molid's into account
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = slice_box(atom,limits) % Basic input arguments
% # atom = slice_box(atom,limits,1) % Inverts the selected sites
%
function atom = slice_box(atom,limits,varargin)

if numel(limits)==3
    limits(4)=limits(1);
    limits(5)=limits(2);
    limits(6)=limits(3);
    limits(1)=0;
    limits(2)=0;
    limits(3)=0;
end

indxlo=find([atom.x]<limits(1));
indxhi=find([atom.x]>limits(4));

indylo=find([atom.y]<limits(2));
indyhi=find([atom.y]>limits(5));

indzlo=find([atom.z]<limits(3));
indzhi=find([atom.z]>limits(6));

ind = unique([indxlo indxhi indylo indyhi indzlo indzhi]);
% molid=unique([atom(intersect([atom.index],ind)).molid]);
% ind = ismember([atom.molid],molid);

if nargin == 3 && cell2mat(varargin(1))==1
    ind_all=1:length([atom.x]);
    ind=setdiff(ind_all,ind);
end

atom(ind)=[];

atom=update_atom(atom);

end

