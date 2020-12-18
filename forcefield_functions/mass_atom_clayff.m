%% mass_atom_clayff.m
% * This function fetches the atom weight from the clayff and interface ff
% * Untested for quite some time
%
%% Version
% 2.081
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = mass_atom_clayff(atom)
% # atom = mass_atom_clayff(atom,'OW')
%
function atom = mass_atom_clayff(atom,varargin)
% Todo... write more general function...

if nargin>1;
    temp_atom=element_atom(atom,'OW','HW',1);
else
    temp_atom=element_atom(atom);
end

clayff_atom=atom;
clayff_param(unique([temp_atom.type]),'spc');
% Check the charge after AssignClayff.m
for i=1:length(temp_atom)
    if strncmpi([temp_atom(i).type],{'Hw'},2);
        ind=strncmpi({'Hw'},[forcefield.clayff.type],2);
    else
        ind=strcmpi([temp_atom(i).type],[forcefield.clayff.type]);
    end
    clayff_atom(i).mass=[forcefield.clayff(ind).mass];
end
%     atom = charge_interface_atom(atom,Box_dim,varargin);

interface_atom=atom;
interface_param(unique([temp_atom.type]),'tip3p');
for i=1:length(temp_atom)
    if strncmpi([temp_atom(i).type],{'Hw'},2);
        ind=strncmpi({'Hw'},[forcefield.interface.type],2);
    else
        ind=strcmp([temp_atom(i).type],[forcefield.interface.type]);
    end
    interface_atom(i).mass=[forcefield.interface(ind).mass];
end

charge_difference=[interface_atom.charge]-[clayff_atom.charge]
if charge_difference == 0;
    [atom.mass]=clayff_atom.mass;
else
    disp('Inconsistencies in the mass assignment between clayff and interface')
    disp('this might still be ok though... just make sure to check the output')
end
