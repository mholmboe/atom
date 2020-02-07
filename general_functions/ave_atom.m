%% ave_atom.m
% * This function calculates the mean of the atom x,y,z coordinates
% * atom is the atom struct
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = ave_atom(atom)
%
function atom = ave_atom(atom)

for i=1:max([atom(:).molid])
    ind=find([atom.molid]==i);
    [atom(ind).ave_x]=deal(mean([atom(ind).x]));
    [atom(ind).ave_y]=deal(mean([atom(ind).y]));
    [atom(ind).ave_z]=deal(mean([atom(ind).z]));
end

