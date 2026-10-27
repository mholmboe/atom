%% check_interface_charge.m
% * This function checks the charge of interface_param atomtypes
% * atom is the atom struct
%
%% Version
% 3.0.1
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = check_interface_charge(atom)
%
%% See also
% # <interface_param.html interface_param>
% # <charge_interface_atom.html charge_interface_atom>
% # <check_interface15_charge.html check_interface15_charge>
% # <charge_atom.html charge_atom>
% # <charge_interface15_atom.html charge_interface15_atom>
%
function atom = check_interface_charge(atom)

nAtoms=size(atom,2);
interface_param(unique([atom.type]),'tip3p');
for i=1:length(atom)
    if strncmpi([atom(i).type],{'Hw'},2);
        ind=strncmpi({'Hw'},[forcefield.interface.type],2);
        atom(i).charge=[forcefield.interface(ind).charge];
    elseif sum(strcmpi([atom(i).type],[forcefield.interface.type])) > 0;
        ind=strcmpi([atom(i).type],[forcefield.interface.type]);
        atom(i).charge=[forcefield.interface(ind).charge];
    else 
        atom(i).charge=0;
    end
end

disp('Total charge')
Total_charge=sum([atom.charge])

assignin('caller','Total_charge',Total_charge);