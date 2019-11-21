%% median_atom.m
% * This function calculates the median position of the atom struct. Should you wrap the atom struct?
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = median_atom(atom)
%
function atom = median_atom(atom)

for i=1:max([atom(:).molid])
    ind=find([atom.molid]==i);
    [atom(ind).med_x]=deal(median([atom(ind).x]));
    [atom(ind).med_y]=deal(median([atom(ind).y]));
    [atom(ind).med_z]=deal(median([atom(ind).z]));
end
