%% interface15_phosphate_atom.m
% * This function tries to assign all atoms according to the interface15 ff
% atom types for hydroxyapatite...
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = interface15_phosphate_atom(atom,Box_dim,1.1,1.8)


function atom = interface15_phosphate_atom(atom,Box_dim,varargin)

if nargin > 2
    rmaxshort=cell2mat(varargin(1));
    rmaxlong=cell2mat(varargin(2));
else
    rmaxshort=1.25; % Maximum allowed H-O bond distance
    rmaxlong=2.25; % Maximum allowed M-O bond distance
end

atom = element_atom(atom);
[atom.element] = atom.type;
[atom.molid]=deal(1);

atom=bond_angle_atom(atom,Box_dim,rmaxshort,rmaxlong,'more'); % Scan all bonds and angles

for i=1:size(atom,2)
    if strncmpi([atom(i).type],'H',1)
        [atom(i).type]=deal({'H'});
        [atom(i).fftype]=deal({'IHOP'});
    end
    if strncmpi([atom(i).type],'P',1)
        [atom(i).type]=deal({'P'});
        [atom(i).fftype]=deal({'IPAP'});
    end
    if strncmpi([atom(i).type],'Ca',2)
        [atom(i).type]=deal({'Ca'});
        [atom(i).fftype]=deal({'ICA_H'});
    end
end

for i=1:size(atom,2)
    Neighbours=sort([atom(i).neigh.type])';
    if numel(Neighbours)>0
        Neighbours=strcat(Neighbours{:});
        if strncmpi([atom(i).type],'O',1)
            if contains(Neighbours,'P')
                [atom(i).type]=deal({'Op'});
                [atom(i).fftype]=deal({'IOAP1'});
            elseif contains(Neighbours,'H')
                [atom(i).type]=deal({'Oh'});
                [atom(i).fftype]=deal({'IOAP2'});
            end
        end
    end
end

atom=bond_angle_atom(atom,Box_dim,rmaxshort,rmaxlong,'more'); % Scan all bonds and angles

atom = check_interface15_charge(atom,'PHOSPHATE');

composition_atom(atom);

