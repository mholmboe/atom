%% tweak_charge_atom.m
% * This function tries to tweak the charge of the atom struct in case of 
% * rounding errors.
% * Please report bugs to michael.holmboe@umu.se
% * Tested 15/4/2017

%% Examples
% * atom=tweak_charge_atom(atom)

function atom=tweak_charge_atom(atom)

nAtoms=size(atom,2);

disp('Total charge before tweaking')
Total_charge=sum([atom.charge])

if abs(round(sum([atom.charge]))-sum([atom.charge])) > 0.000001 && abs(round(sum([atom.charge]))-sum([atom.charge])) < 0.5
    disp('Tweaking the charge')
    qtot=sum([atom.charge]);
    delta_q=sum([atom.charge])-round(sum([atom.charge]));
    charge=num2cell([atom.charge]-delta_q/nAtoms); [atom.charge]=deal(charge{:});
    Total_charge=sum([atom.charge])
end

assignin('caller','Total_charge',Total_charge);

