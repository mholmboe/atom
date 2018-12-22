%% tweak_charge_atom.m
% * This function tries to tweak the charge of the atom struct in case of
% rounding errors.
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=tweak_charge_atom(atom)
%
function atom=tweak_charge_atom(atom)

nAtoms=size(atom,2);

disp('Total charge before tweaking')
Total_charge=sum([atom.charge])

for i=1:size(atom,2)
    atom(i).charge=round(atom(i).charge,5);
end

if abs(round(sum([atom.charge]))-sum([atom.charge])) > 0.00001 && abs(round(sum([atom.charge]))-sum([atom.charge])) < 0.5
    disp('Tweaking the charge')
    qtot=sum([atom.charge]);
    delta_q=sum([atom.charge])-round(sum([atom.charge]));
    if abs(delta_q) < 0.05
        [atom(1).charge]=atom(1).charge-delta_q;
    else
        charge=num2cell([atom.charge]-delta_q/nAtoms); [atom.charge]=deal(charge{:});
        pause(5)
    end
    Total_charge=sum([atom.charge])
end

assignin('caller','Total_charge',Total_charge);

