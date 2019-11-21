%% check_clayff_2004_charge.m 
% * This function checks the charge of clayff 2004 atomtypes
% * atom is the atom struct
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = check_clayff_charge(atom)
%
function atom = check_clayff_2004_charge(atom)

nAtoms=size(atom,2);
clayff_2004_param(unique([atom.type]),'SPC/E');
% Check the charge after clayff_atom.m
for i=1:length(atom)
    if strncmpi([atom(i).type],{'Hw'},2)
        ind=strncmpi({'Hw'},[forcefield.clayff.type],2);
        atom(i).charge=[forcefield.clayff(ind).charge];
    elseif sum(strcmpi([atom(i).type],[forcefield.clayff.type])) > 0
        ind=strcmpi([atom(i).type],[forcefield.clayff.type]);
        atom(i).charge=[forcefield.clayff(ind).charge];
    else 
        atom(i).charge=0;
    end
end
disp('Total charge')
Total_charge=sum([atom.charge])

assignin('caller','Total_charge',Total_charge);
