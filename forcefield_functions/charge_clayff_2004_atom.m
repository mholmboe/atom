%% charge_clayff_2004_atom.m
% * This function tries to smear out the charge at isomorphic substitutions sites according to clayff 
% * atom is the atom struct
% * Box_dim is the box dimension vector
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples with some fixed charges, the rest is smeared over the O's
% # Total_charge = charge_clayff_atom(atom,Box_dim,{'ao' 'mgo' 'st' 'ho'},[1.575 1.36 2.1 0.425])

function atom = charge_clayff_2004_atom(atom,Box_dim,varargin)
nAtoms=size(atom,2);
[atom.charge]=deal(0);

% Atom_label=sort(unique([atom.type]));
% clayff_2004_param(sort(Atom_label),watermodel);
% no_O_label=Atom_label(~strncmp(Atom_label,'O',1));
% no_O_ind=ismember(Atom_label,no_O_label);
% atom = charge_clayff_2004_atom(atom,Box_dim,Atom_label(no_O_ind),Charge(no_O_ind));

if nargin>2
    Atom_label=varargin{1}(:);
    Charge=cell2mat(varargin(2));
    [Atom_label,sort_ind]=sort(Atom_label);
    Charge=Charge(sort_ind);
    Met_ind=zeros(1,nAtoms);
    for i=1:length(Charge)
        ind=strcmpi([atom.type],Atom_label(i));
        [atom(ind).charge]=deal(Charge(i));
        Met_ind=Met_ind+ind;
    end
    Met_ind=find(Met_ind);
    Ox_ind=setdiff(1:nAtoms,Met_ind);
    
    atom=bond_angle_atom(atom,Box_dim,1.25,2.1,'more');

    for i=1:length(Ox_ind)
        bond_ind=setdiff(reshape(atom(Ox_ind(i)).bond.index,[],1),Ox_ind(i));
        Zsum=0;
        if ~isempty(bond_ind)
            if bond_ind(1)>0
                for j=1:length(bond_ind)
                    if strncmpi([atom(bond_ind(j)).type],'st',2)
                        Z=4;
                    elseif strncmpi([atom(bond_ind(j)).type],'ao',2)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'feo',2)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'mg',2)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'mn',2)
                        Z=4;
                    elseif strncmpi([atom(bond_ind(j)).type],'cao',3)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'ho',1)
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
    clayff_2004_param(unique([atom.type]),'SPC/E');
    for i=1:length(atom)
        if strncmpi([atom(i).type],{'Hw'},2)
            ind=strncmpi({'Hw'},[forcefield.clayff.type],2);
        else
            ind=strcmpi([atom(i).type],[forcefield.clayff.type]);
        end
        atom(i).charge=[forcefield.clayff(ind).charge];
    end
end

disp('Total charge')
Total_charge=sum([atom.charge])

assignin('caller','Total_charge',Total_charge);


