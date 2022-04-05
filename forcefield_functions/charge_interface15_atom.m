%% charge_interface15_charge.m
% * This function tries to smear out the charge at isomorphic substitutions sites sites according to the interface ff
% * atom is the atom struct
% * Box_dim is the box dimension vector
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples with some fixed charges, the rest is smeared over the O's
% # Total_charge = charge_interface15_atom(atom,Box_dim,{'AY1' 'MY1' 'SY1' 'HOY'},[1.45 1.1 1.1 0.4])

function atom = charge_interface15_atom(atom,Box_dim,varargin)

nAtoms=size(atom,2);
[atom.charge]=deal(0);

% Atom_label=sort(unique([atom.type]));
% clayff_param(sort(Atom_label),watermodel);
% no_O_label=Atom_label(~strncmp(Atom_label,'O',1));
% no_O_ind=ismember(Atom_label,no_O_label);
% atom = charge_clayff_atom(atom,Box_dim,Atom_label(no_O_ind),Charge(no_O_ind));

if nargin>2
    Atom_label=varargin{1}(:);
    Charge=cell2mat(varargin(2));
    [Atom_label,sort_ind]=sort(Atom_label);
    Charge=Charge(sort_ind);
    Met_ind=zeros(1,nAtoms);
    for i=1:length(Charge)
        ind=strcmpi([atom.fftype],Atom_label(i));
        [atom(ind).charge]=deal(Charge(i));
        Met_ind=Met_ind+ind;
    end
    Met_ind=find(Met_ind);
    Ox_ind=setdiff(1:nAtoms,Met_ind);
    
    
    atom=bond_angle_atom(atom,Box_dim,1.25,2.25,'more');
    for i=1:length(Ox_ind)
        bond_ind=setdiff(reshape(atom(Ox_ind(i)).bond.index,[],1),Ox_ind(i));
        Zsum=0;
        if ~isempty(bond_ind)
            if bond_ind(1)>0
                for j=1:length(bond_ind)
                    if strncmpi([atom(bond_ind(j)).type],'Si',2)
                        Z=4;
                    elseif strncmpi([atom(bond_ind(j)).type],'Al',2)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'Feo',2)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'Mg',2)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'H',1)
                        Z=1;
                    else
                        Z=0;
                    end
                    Zp=atom(bond_ind(j)).charge;
                    CN=size(atom(bond_ind(j)).bond.index,1);
                    Zsum=Zsum+(Z-Zp)/CN;
                end
            end
        end
        atom(Ox_ind(i)).charge= -2.00 + Zsum;
    end
else
    interface15_param(unique([atom.fftype]),'tip3p');
    for i=1:length(atom)
        if strncmpi([atom(i).fftype],{'Hw'},2);
            ind=strncmpi({'Hw'},[forcefield.interface.type],2);
        else
            ind=strcmpi([atom(i).fftype],[forcefield.interface.type]);
        end
        atom(i).charge=[forcefield.interface(ind).charge];
    end
end
Total_charge=sum([atom.charge])
assignin('caller','Total_charge',Total_charge);


