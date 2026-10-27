%% median_atom.m
% * This function calculates the median position of the atom struct. Should you wrap the atom struct?
%
%% Version
% 3.0.2
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = median_atom(atom)
%
%% See also
% # <ave_atom.html ave_atom>
% # <COM_atom.html COM_atom>
% # <center_atom.html center_atom>
% # <place_atom.html place_atom>
%
function atom = median_atom(atom)

for i=1:max([atom(:).molid])
    ind=find([atom.molid]==i);
    [atom(ind).med_x]=deal(median([atom(ind).x]));
    [atom(ind).med_y]=deal(median([atom(ind).y]));
    [atom(ind).med_z]=deal(median([atom(ind).z]));
end
