%% interface15_silica_atom.m
% * This function tries to assign all atoms according to the interface15 ff
% atom types for Silica made up by me...
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = interface15_silica_atom(atom,Box_dim,1.1,1.8)


function atom = interface15_silica_atom(atom,Box_dim,varargin)

if nargin > 2
    rmin=cell2mat(varargin(1));
    rlarge=cell2mat(varargin(2));
else
    rmin=1.1; % Maximum allowed H-O bond distance
    rlarge=2.1; % Maximum allowed M-O bond distance
end

atom = element_atom(atom);
[atom.element] = atom.type;
[atom.molid]=deal(1);

atom=bond_angle_atom(atom,Box_dim,rmin,rlarge,'more'); % Scan all bonds and angles

for i=1:size(atom,2)
    if strncmpi([atom(i).type],'H',1)
        [atom(i).type]=deal({'H'});
        [atom(i).fftype]=deal({'HOY'});
    end
end

for i=1:size(atom,2)
    i;
    Neighbours=sort([atom(i).neigh.type])';
    Neighbours=strcat(Neighbours{:});
    if strncmpi([atom(i).type],'O',1)
        if contains(Neighbours,'SiSi')
            [atom(i).type]=deal({'Ob'});
            [atom(i).fftype]=deal({'OC23'});
        elseif contains(Neighbours,'HSi')
            [atom(i).type]=deal({'Osih'});
            [atom(i).fftype]=deal({'OC24'});
        elseif strncmpi(Neighbours,'Si',2)
            [atom(i).type]=deal({'Osi'});
            [atom(i).fftype]=deal({'OC25'});
        end
    end
end

atom=bond_angle_atom(atom,Box_dim,rmin,rlarge,'more'); % Scan all bonds and angles

for i=1:size(atom,2)
    if strncmpi([atom(i).type],'Si',1)
        Neighbours=sort([atom(i).neigh.type])';
        Neighbours=strcat(Neighbours{:});
        [atom(i).type]=deal({'Si'});
        [atom(i).fftype]=deal({'SC4'});
        if contains(Neighbours,'Osi')
            [atom(i).type]=deal({'Si'});
            [atom(i).fftype]=deal({'SC5'});
        end
        Neigh_ind=[atom(i).neigh.index];
        Second_neighbours=0;
        for nn=1:numel(Neigh_ind)
            Second_neighbours=Second_neighbours+length([atom(Neigh_ind(nn)).neigh.index]);
        end
        %         if Second_neighbours==8
        %             atom(i).fftype={'SC4'};
        %             atom(i).type={'Si'};
        %         elseif Second_neighbours==7
        %             atom(i).fftype={'SC5'};
        %             atom(i).type={'Si'};
        %         else
        if Second_neighbours < 7
            atom(i).fftype={'SC5'};
            atom(i).type={'Si'};
            disp('Si has more than one undercoordinated oxygen neighbours')
            Second_neighbours
            i
        end
    end
end

atom = check_interface15_charge(atom,'SILICA');

composition_atom(atom);

