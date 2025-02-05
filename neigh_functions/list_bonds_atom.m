%% list_bonds_angles_atom.m
% * This function tries to find all bonds and angles between the  atomtypes
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # [Bond_data,Angle_data] = list_bonds_atom(atom,Box_dim)
% # [Bond_data,Angle_data] = list_bonds_atom(atom,[],Bond_index)
% # [Bond_data,Angle_data] = list_bonds_atom(atom,[],Bond_index,Angle_index)
% # [Bond_data,Angle_data,Dihedral_data] = list_bonds_atom(atom,[],Bond_index,Angle_index,Dihedral_index)

function [Bond_data,Angle_data] = list_bonds_atom(atom,Box_dim,varargin)

Bond_data=[];Angle_data=[];
if numel(Box_dim)>0
    if nargin > 2
        rmaxshort=varargin{1};
        rmaxlong=varargin{2};
    else
        rmaxshort=1.25;
        rmaxlong=2.45;
    end
    atom=bond_angle_atom(atom,Box_dim,rmaxshort,rmaxlong);
    Bond_data = map_bonded(atom,Bond_index(:,1:3));
    Angle_data = map_bonded(atom,Angle_index(:,1:4));
    % Dihedral_data = map_bonded(atom,Dihedral_index(:,1:5));
end

if nargin>2
    Bond_index=varargin{1};
    Bond_data = map_bonded(atom,Bond_index(:,1:3));
end
if nargin>3
    Angle_index=varargin{2};
    Angle_data = map_bonded(atom,Angle_index(:,1:4));
end

end

function result = map_bonded(atom,index)

% Initialize a containers.Map to store bond distances for each atom pair
Map = containers.Map('KeyType', 'char', 'ValueType', 'any');

Atom_labels=unique([atom.type]);
A_labels=[Atom_labels(strncmp(Atom_labels,'O',1)) Atom_labels(strncmp(Atom_labels,'Fs',2))];
M_labels=Atom_labels(~ismember(Atom_labels,A_labels));

if size(index,2)==3
    % Iterate over each row in the bonds matrix
    for i = 1:size(index, 1)
        % Get atom indices
        atom1 = index(i,1);
        atom2 = index(i,2);

        % Get corresponding atom types
        atomType1 = string(atom(atom1).type);
        atomType2 = string(atom(atom2).type);

        % Sort atom types alphabetically to avoid duplicate pairs
        if sum(ismember(atomType1,M_labels))>0
            pairKey = atomType1 + atomType2;
        else
            pairKey = atomType2 + atomType1;
        end

        % Get bond distance
        Val = index(i, 3);

        % Store bond distance in the map
        if isKey(Map, pairKey)
            Map(pairKey) = [Map(pairKey), Val];
        else
            Map(pairKey) = Val;
        end
    end

elseif size(index,2)==4

    % Iterate over each row in the bonds matrix
    for i = 1:size(index, 1)
        % Get atom indices
        atom1 = index(i,1);
        atom2 = index(i,2);
        atom3 = index(i,3);

        % Get corresponding atom types
        atomType1 = string(atom(atom1).type);
        atomType2 = string(atom(atom2).type);
        atomType3 = string(atom(atom3).type);

        % Sort atom types alphabetically to avoid duplicate pairs
        if atomType1 < atomType3
            pairKey = atomType1 + atomType2 + atomType3;
        else
            pairKey = atomType3 + atomType2 + atomType1;
        end

        % Get bond distance
        Val = index(i, 4);

        % Store bond distance in the map
        if isKey(Map, pairKey)
            Map(pairKey) = [Map(pairKey), Val];
        else
            Map(pairKey) = Val;
        end
    end

elseif size(index,2)==5

    % % Iterate over each row in the bonds matrix
    % for i = 1:size(index, 1)
    %     % Get atom indices
    %     atom1 = index(i,1);
    %     atom2 = index(i,2);
    %     atom3 = index(i,3);
    %     atom3 = index(i,4);
    % 
    %     % Get corresponding atom types
    %     atomType1 = string(atom(atom1).type);
    %     atomType2 = string(atom(atom2).type);
    %     atomType3 = string(atom(atom3).type);
    %     atomType4 = string(atom(atom4).type);
    % 
    %     % Sort atom types alphabetically to avoid duplicate pairs
    %     if atomType1 < atomType4
    %         pairKey = atomType1 + atomType2 + atomType3 + atomType4;
    %     else
    %         pairKey = atomType4 + atomType3 + atomType2 + atomType1;
    %     end
    % 
    %     % Get bond distance
    %     Val = index(i, 5);
    % 
    %     % Store bond distance in the map
    %     if isKey(Map, pairKey)
    %         Map(pairKey) = [Map(pairKey), Val];
    %     else
    %         Map(pairKey) = Val;
    %     end
    % end

end

% Prepare output matrix to store the results
result = cell(length(Map), 2);
i = 1;

% Calculate average bond distance for each pair and store in result
for key = keys(Map)
    Values = Map(key{1});
    averageValue = mean(Values);
    stdValue = std(Values);
    result{i, 1} = key{1};
    result{i, 2} = averageValue;
    result{i, 3} = stdValue;
    i = i + 1;
end

end