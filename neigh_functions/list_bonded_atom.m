%% list_bonded_atom.m
% * This function tries to find all bonds, angles or dihedral between the  atomtypes
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # [Bonded_list] = list_bonds_atom(atom,Bond_index)

function Bond_order = list_bonded_atom(atom,Bond_index)


if size(Bond_index,2)==3
    %% Bonds

    N = size(Bond_index, 1);
    bondDistances = Bond_index(:,3);

    bondTypes = cell(N, 2);
    for k = 1:N
        idx1 = Bond_index(k, 1);
        idx2 = Bond_index(k, 2);
        atype1 = char(atom(idx1).type);
        atype2 = char(atom(idx2).type);
        sortedPair = sort({atype1, atype2});
        bondTypes(k,:) = sortedPair;
    end

    % N = size(bondTypes, 1);
    pairsCombined = cell(N, 1);

    for i = 1:N
        a = bondTypes{i,1};
        b = bondTypes{i,2};
        pairsCombined{i} = [a, '_', b];
    end

    [uniquePairsCombined, ~, ~] = unique(pairsCombined);

    avgDistances = zeros(numel(uniquePairsCombined),1);
    stdDistances = zeros(numel(uniquePairsCombined),1);
    for i = 1:numel(uniquePairsCombined)
        matches = strcmp(pairsCombined, uniquePairsCombined{i});
        avgDistances(i) = mean(bondDistances(matches));
        stdDistances(i) = 100*(std(bondDistances(matches))./mean(bondDistances(matches)));
    end

    uniquePairsCell = cellfun(@(x) strsplit(x,'_'), uniquePairsCombined, 'UniformOutput', false);
    uniquePairs = vertcat(uniquePairsCell{:});
    Bond_order = [uniquePairs, num2cell(avgDistances), num2cell(stdDistances)];

end

if size(Bond_index,2)==4
    Angle_index=Bond_index;

    %% Angles

    angleTypes = cell(size(Angle_index,1), 3);
    angleAngles = Angle_index(:,4);

    for i = 1:size(Angle_index,1)
        idx1 = Angle_index(i,1);
        idx2 = Angle_index(i,2);
        idx3 = Angle_index(i,3);
        type1 = char(atom(idx1).type);
        type2 = char(atom(idx2).type);
        type3 = char(atom(idx3).type);
        sortedEnds = sort({type1, type3});
        angleTypes(i,:) = {sortedEnds{1}, type2, sortedEnds{2}};
        % angleTypes(i,:) = {type1, type2, type3};
    end

    N = size(angleTypes, 1);
    tripletsCombined = cell(N, 1);

    for i = 1:N
        a = angleTypes{i,1};
        b = angleTypes{i,2};
        c = angleTypes{i,3};
        tripletsCombined{i} = [a, '_', b,'_', c];
    end

    [uniqueTripletsCombined, ~, ~] = unique(tripletsCombined);

    avgAngles = zeros(numel(uniqueTripletsCombined),1);
    stdAngles = zeros(numel(uniqueTripletsCombined),1);
    for i = 1:numel(uniqueTripletsCombined)
        matches = strcmp(tripletsCombined, uniqueTripletsCombined{i});
        avgAngles(i) = mean(angleAngles(matches));
        stdAngles(i) = 100*(std(angleAngles(matches))./mean(angleAngles(matches)));
    end

    uniqueTripletsCell = cellfun(@(x) strsplit(x,'_'), uniqueTripletsCombined, 'UniformOutput', false);
    uniqueTriplets = vertcat(uniqueTripletsCell{:});
    Angle_order = [uniqueTriplets, num2cell(avgAngles), num2cell(stdAngles)];

    Bond_order=Angle_order;

end

% if size(Bond_index,2)==5
%     Dihedral_index=Bond_index;
%     %% Dihedrals
% 
%     % Preallocate a cell array for storing dihedral type quadruples
%     dihedralTypes = cell(size(Dihedral_index,1), 4);
% 
%     % Extract dihedral types based on atom indices
%     for i = 1:size(Dihedral_index,1)
%         % Get the atom indices
%         idx1 = Dihedral_index(i,1);
%         idx2 = Dihedral_index(i,2);
%         idx3 = Dihedral_index(i,3);
%         idx4 = Dihedral_index(i,4);
% 
%         % Get the corresponding atom types
%         type1 = char(atom(idx1).type);
%         type2 = char(atom(idx2).type);
%         type3 = char(atom(idx3).type);
%         type4 = char(atom(idx4).type);
% 
%         % Check the alphabetical order of type1 and type4
%         firstind=find(strcmp(sort([type1, type4]),type1));
%         if firstind==1
%             % If type1 is alphabetically first or equal, keep order as is:
%             % (type1, type2, type3, type4)
%             quadruple = {type1, type2, type3, type4};
%         else
%             % If type4 should come first alphabetically, reverse the order:
%             % (type4, type3, type2, type1)
%             quadruple = {type4, type3, type2, type1};
%         end
% 
%         % Store the standardized quadruple
%         dihedralTypes(i,:) = quadruple;
%     end
% 
%     % Find all unique dihedral quadruples
%     uniqueQuadruples = unique(dihedralTypes, 'rows');
% 
%     % Initialize the Dihedral_order cell array with five columns
%     % (four atom types + average dihedral angle)
%     Dihedral_order = cell(size(uniqueQuadruples,1),5);
% 
%     % Compute the average dihedral angle for each unique quadruple
%     for j = 1:size(uniqueQuadruples,1)
%         % Find the rows in dihedralTypes that match this unique quadruple
%         match_rows = ismember(dihedralTypes, uniqueQuadruples(j,:),'rows');
% 
%         % Average the dihedral angles for these matched rows
%         avg_dihedral = mean(Dihedral_index(match_rows,5));
% 
%         % Store the result: four atom types and the average dihedral
%         Dihedral_order(j,:) = [uniqueQuadruples(j,:), {avg_dihedral}];
%     end
%     Bond_order=Dihedral_order;
% end

end