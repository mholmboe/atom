function atom = charge_opls_go_atom(atom,Box_dim,varargin)
%% This is a tailor-made function and likely not relvant for you!!!
%%
%% charge_opls_go_atom.m - This function tries to smear out the charge at around -OH and epoxides in GO
%% Based on the corresponding clayff atom function
%% Use like this to set new charges
%% Total_charge = charge_clayff_atom(atom,Box_dim,{'Al' 'Mgo' 'Si' 'H'},[1.575 1.36 2.1 0.425])

nAtoms=size(atom,2);
[atom.charge]=deal(0);

% Atom_label=sort(unique([atom.type]));
% clayff_param(sort(Atom_label),watermodel);
% no_O_label=Atom_label(~strncmp(Atom_label,'O',1));
% no_O_ind=ismember(Atom_label,no_O_label);
% atom = charge_clayff_atom(atom,Box_dim,Atom_label(no_O_ind),Charge(no_O_ind));

atom_C=atom;
%if nargin>2;
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

atom=bond_angle_atom(atom,Box_dim,1.05,1.6,'more');
for i=1:length(Ox_ind)
    bond_ind=setdiff(reshape(atom(Ox_ind(i)).bond.index,[],1),Ox_ind(i));
    Zsum=0;
    if ~isempty(bond_ind)
        if bond_ind(1)>0
            for j=1:length(bond_ind)
                if strncmpi([atom(bond_ind(j)).type],'H',2)
                    Z=0;
                    CN=1;
                    Zp=atom(bond_ind(j)).charge;
                elseif strncmpi([atom(bond_ind(j)).type],'Oe',2)
                    Z=0;
                    CN=2;
                    charge_ind=setdiff(reshape(atom(bond_ind(j)).bond.index,[],1),bond_ind(j));
                    Zp=atom(bond_ind(j)).charge;%+sum([atom(charge_ind).charge]);
                elseif strncmpi([atom(bond_ind(j)).type],'Oh',2)
                    Z=0;
                    CN=1;
                    charge_ind=setdiff(reshape(atom(bond_ind(j)).bond.index,[],1),bond_ind(j));
                    Zp=atom(bond_ind(j)).charge+sum([atom(charge_ind).charge]);
                elseif strncmpi([atom(bond_ind(j)).type],'C',1)
                    Z=0;
                    CN=4;
                    Zp=atom(bond_ind(j)).charge;
                else
                    Z=0;
                end
                %                                          Zp=atom(bond_ind(j)).charge;
                %                                          CN=size(atom(bond_ind(j)).bond.index,1);
                Zsum=Zsum+(Z-Zp)/CN;
            end
        end
    end
    atom_C(Ox_ind(i)).charge= Zsum;
end
% else
%     clayff_param(unique([atom.type]),'SPC/E');
%     %% Check the charge after AssignClayff.m
%     for i=1:length(atom)
%         if strncmpi([atom(i).type],{'Hw'},2);
%             ind=strncmpi({'Hw'},[forcefield.clayff.type],2);
%         else
%             ind=strcmpi([atom(i).type],[forcefield.clayff.type]);
%         end
%         atom(i).charge=[forcefield.clayff(ind).charge];
%     end
%end

atom(Ox_ind)=atom_C(Ox_ind);

for i=1:size(atom,2)
    if strcmp(atom(i).type,'Cen')
        for j=1:length([atom(i).neigh.index])
            %             i
            %             j
            if sum(strcmp(atom(i).neigh.type(j),{'Ce' 'Coh'})) > 0
                %                 i
                %                 j
                [atom(i).charge]=[atom(i).charge]+0.02;
                [atom([atom(i).neigh.index(j)]).charge]=[atom([atom(i).neigh.index(j)]).charge]-0.02;
            end
        end
    end
end


nAtoms=size(atom,2);
disp('Total charge without tweaking')
Total_charge=sum([atom.charge])

if nargin>4;
    disp('Tweaking the charges of all atoms with charge >= +0.1')
    qtot=sum([atom.charge]);
    delta_q=sum([atom.charge])-0;%round(sum([atom.charge]));
    ind_high_charge=find([atom([atom.charge]>+0.1).index]);
    nhigh_charge=length(ind_high_charge);
    charge_C=num2cell([atom(ind_high_charge).charge]-delta_q/nhigh_charge); [atom(ind_high_charge).charge]=deal(charge_C{:});
    
    disp('Total charge after tweaking')
    Total_charge=sum([atom.charge])
end


% atom=tweak_charge_atom(atom)

% assignin('caller','atom_C',atom_C);
assignin('caller','Total_charge',Total_charge);



