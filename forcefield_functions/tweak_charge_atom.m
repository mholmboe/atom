%% tweak_charge_atom.m
% * This function tries to tweak the charge of the atom struct in case of
% rounding errors.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=tweak_charge_atom(atom)
% # atom=tweak_charge_atom(atom,{'Si'}) % Will only tweak the charges for the atomtype 'Si', and make the system neutral
%
function atom=tweak_charge_atom(atom,varargin)

nAtoms=size(atom,2);

disp('Total charge before tweaking')
Total_charge=sum([atom.charge])

for i=1:size(atom,2)
    atom(i).charge=round(atom(i).charge,6);
end

if nargin>1
    Atom_labels=varargin{1};
    ind=ismember([atom.type],Atom_labels);
    nAtoms=size(atom(ind),2);
    disp('Tweaking the charge for the atomtype')
    Atom_labels
    qtot=sum([atom.charge]);
    delta_q=sum([atom.charge]);%-round(sum([atom.charge]));
    charge=num2cell([atom(ind).charge]-delta_q/nAtoms);
    [atom(ind).charge]=deal(charge{:});
    Total_charge=sum([atom.charge])
else
    if abs(round(sum([atom.charge]))-sum([atom.charge])) > 0.00001 && abs(round(sum([atom.charge]))-sum([atom.charge])) < 0.5
        disp('Tweaking the charge')
        qtot=sum([atom.charge]);
        delta_q=sum([atom.charge])-round(sum([atom.charge]));
        if abs(delta_q) < 0.05
            [atom(1).charge]=atom(1).charge-delta_q;
        else
            charge=num2cell([atom.charge]-delta_q/nAtoms);
            [atom.charge]=deal(charge{:});
        end
        Total_charge=sum([atom.charge])
    end
end

for i=1:size(atom,2)
    atom(i).charge=round(atom(i).charge,6);
end

final_delta_q=sum([atom.charge])-round(sum([atom.charge]));

if abs(final_delta_q) < 0.05
    [atom(1).charge]=atom(1).charge-final_delta_q;
end

round(unique([atom.charge]),8)

disp('Total charge after tweaking')
Total_charge=sum([atom.charge])

assignin('caller','Total_charge',Total_charge);

end


