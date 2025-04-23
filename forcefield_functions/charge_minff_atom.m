%% charge_minff_atom.m
% * This function tries to smear out the charge at isomorphic substitutions
% * sites according to minff
%
% * atom is the atom struct
% * Box_dim is the box dimension vector
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples with some fixed charges, the rest is smeared over the O's
% * atom = charge_minff_atom(atom,Box_dim,{'Al' 'Mgo' 'Si' 'H'},[1.575 1.36 2.1 0.425])

function atom = charge_minff_atom(atom,Box_dim,varargin)
nAtoms=size(atom,2);
[atom.charge]=deal(0);

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
    Ox_ind=find(strncmpi([atom.type],'O',1)); % setdiff(1:nAtoms,Met_ind);
    Fs_ind=find(strncmpi([atom.type],'Fs',2)); % setdiff(1:nAtoms,Met_ind);

    % atom = bond_angle_atom(atom,Box_dim,1.25,2.45);%,'more');
    atom = bond_atom(atom,Box_dim,2.45,0.6);

    if numel(Ox_ind) >= size(Bond_index,1)
        atom = bond_atom(atom,Box_dim,2.45,0.6);
        size(Bond_index,1)
    end

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
                    elseif strncmpi([atom(bond_ind(j)).type],'Fe2',3)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'Feo2',4)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'Fet2',4)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'Fee2',4)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'Feo3',4)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'Fet3',4)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'Fee3',4)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'Fe',2)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'F',1) % Fs
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'Ti',2)
                        Z=4;
                    elseif strncmpi([atom(bond_ind(j)).type],'Li',2)
                        Z=1;
                    elseif strncmpi([atom(bond_ind(j)).type],'Mg',2)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'Ca',2)
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
        atom(Ox_ind(i)).charge = -2.00 + Zsum;
    end

    for i=1:length(Fs_ind)
        bond_ind=setdiff(reshape(atom(Fs_ind(i)).bond.index,[],1),Fs_ind(i));
        Zsum=0;
        if ~isempty(bond_ind)
            if bond_ind(1)>0
                for j=1:length(bond_ind)
                    if strncmpi([atom(bond_ind(j)).type],'Si',2)
                        Z=4;
                    elseif strncmpi([atom(bond_ind(j)).type],'Al',2)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'Fe2',3)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'Feo2',4)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'Fet2',4)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'Fee2',4)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'Feo3',4)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'Fet3',4)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'Fee3',4)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'Fe',2)
                        Z=3;
                    elseif strncmpi([atom(bond_ind(j)).type],'Ti',2)
                        Z=4;
                    elseif strncmpi([atom(bond_ind(j)).type],'Li',2)
                        Z=1;
                    elseif strncmpi([atom(bond_ind(j)).type],'Mg',2)
                        Z=2;
                    elseif strncmpi([atom(bond_ind(j)).type],'Ca',2)
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
        atom(Fs_ind(i)).charge = -1.00 + Zsum;
    end

    [Atom_label,at_ind]=unique([atom.type],'stable');
    %     Atom_label
    Charge=unique([atom(at_ind).charge],'stable');
end

ind_Hw=strncmpi([atom.type],'Hw',2);
ind_Ow=strncmpi([atom.type],'Ow',2);
[atom(ind_Hw).charge]=deal(0.447585);
[atom(ind_Ow).charge]=deal(-0.89517);

% atom = composition_atom(atom);

disp('Total charge')
Total_charge=sum([atom.charge])
if round2dec(Total_charge)~=sum(Total_charge)
    disp('Run tweak_charge_atom() to get an integer charge of the struct')
end

assignin('caller','Atom_labels',Atom_label);
assignin('caller','Charges',Charge);
assignin('caller','Total_charge',Total_charge);



