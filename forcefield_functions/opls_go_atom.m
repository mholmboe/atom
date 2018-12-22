%% charge_opls_go_atom.m
% * This function tries to smear out the charge at around -OH and epoxides in GO
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom = opls_go_atom(atom,Box_dim,1.1,1.8)


function atom = opls_go_atom(atom,Box_dim,rmin,rlarge)

% atom=atom_full;
% rmin=1.01
% rlarge=1.8

atom = element_atom(atom);
[atom.element] = atom.type;
[atom.molid]=deal(1);
atom=bond_angle_atom(atom,Box_dim,rmin,rlarge,'more'); % Scan all bonds and angles

for i=1:size(atom,2)
    i;
    Neighbours=sort([atom(i).neigh.type])';
    Neighbours=strcat(Neighbours{:});
    if strncmpi([atom(i).type],'O',1)
        if strncmpi(Neighbours,'CH',2)
            %       if strncmpi([atom_full(i).type],'Oh',2)
            [atom(i).type]=deal({'Oh'});
            [atom(i).fftype]=deal({'opls_154'});
        end
        if strncmpi(Neighbours,'CC',2)
            %       if strncmpi([atom_full(i).type],'Oe',2)
            [atom(i).type]=deal({'Oe'});
            [atom(i).fftype]=deal({'opls_180'});
        end
    elseif strncmpi([atom(i).type],'H',1)
        [atom(i).type]=deal({'H'})
        [atom(i).fftype]=deal({'opls_155'});
    end
end

atom=bond_angle_atom(atom,Box_dim,1.01,1.6,'more'); % Scan all bonds and angles

for i=1:size(atom,2)
    Neighbours=sort([atom(i).neigh.type])';
    Neighbours=strcat(Neighbours{:});
    if strncmpi([atom(i).type],'C',1)
        if numel([atom(i).neigh.index]) == 3
            [atom(i).type]=deal({'Cen'});
            [atom(i).fftype]=deal({'opls_141'});
        end
        
        if strncmpi(Neighbours,'CCCOh',5)
            [atom(i).type]=deal({'Coh'});
            [atom(i).fftype]=deal({'opls_159'});
        end
        
        if strncmpi(Neighbours,'CCCOe',5)
            [atom(i).type]=deal({'Ce'});
            [atom(i).fftype]=deal({'opls_183'});
        end
    end
end

atom=bond_angle_atom(atom,Box_dim,1.01,1.8,'more'); % Scan all bonds and angles

