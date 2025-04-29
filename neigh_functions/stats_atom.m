%% stats_atom.m Generate statistics about atom types, coordination, and charges in the structure.
% * This function analyzes atom types, their coordination environment, charges,
% coordination numbers, bond distances, and angles, and outputs a formatted report.
% The report can be written to a log file and/or returned as a string.
%
% Inputs:
%   atom - Structure array containing atom data with fields:
%          .type - Atom type string
%          .neigh - Array of neighbor indices
%          .charge - Atom charge value
%          .bond - Bonds data (if available)
%          .angles - Angles data (if available)
%          .mass - Atom mass (if available)
%   total_charge - The total charge of the system
%   ffname - The name of the forcefield used (e.g., 'minff', 'clayff')
%   Box_dim - Optional box dimensions array. Can be a 3-element array [Lx, Ly, Lz] for
%            orthogonal boxes or a 6/9-element array for triclinic boxes.
%   Cell - Optional cell parameters array [a, b, c, alpha, beta, gamma]
%   log_file - Optional path to a log file. If provided, statistics will be
%             written to this file.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
%   result =  stats_atom(atom, Box_dim)
%   rresult =  stats_atom(atom, Box_dim, log_file)
%
function result =  stats_atom(atom,Box_dim,log_file)

% Box_dim to Cell
if size(Box_dim,2)==6
    Cell=Box_dim;
    Box_dim=Cell2Box_dim(Cell);
else
    Cell=Box_dim2Cell(Box_dim);
end

if nargin < 3
    log_file = 'output.log';
end

% Conversion factor for density calculation: amu/Å³ to g/cm³
% 1 amu = 1.66053886e-24 g and 1 Å³ = 1e-24 cm³
AMU_TO_G_PER_CM3 = 1.66053886;

% Initialize statistics storage
all_neighbor_str_list = {};
atom_type_counts = containers.Map; %('KeyType', 'char', 'ValueType', 'any');
coord_nums = containers.Map('KeyType', 'char', 'ValueType', 'any');
bond_dists = containers.Map('KeyType', 'char', 'ValueType', 'any');
angles = containers.Map('KeyType', 'char', 'ValueType', 'any');
h_involved_angles = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Additional statistics storage for atom type pairs and triplets
bond_type_pairs = containers.Map('KeyType', 'char', 'ValueType', 'any');
angle_type_triplets = containers.Map('KeyType', 'char', 'ValueType', 'any');

%
atom=bond_atom(atom,Box_dim,2.45);

% Gather coordination information
for i = 1:length(atom)
    atom_type = '';
    if isfield(atom, 'type')
        atom_type = string(atom(i).type);
    else
        atom_type = "X";
    end

    % Count atom types
    if isKey(atom_type_counts, atom_type)
        atom_type_counts(atom_type) = atom_type_counts(atom_type) + 1;
    else
        atom_type_counts(atom_type) = 1;
    end

    % Get neighbor atom types for coordination environment
    neighbor_indices = [];
    if isfield(atom, 'neigh')
        neighbor_indices = atom(i).neigh.index;
    end

    neighbor_types = string(atom(i).neigh.type)';

    % Collect coordination number
    cn = length(neighbor_indices);
    if isKey(coord_nums, atom_type)
        coord_nums(atom_type) = [coord_nums(atom_type), cn];
    else
        coord_nums(atom_type) = cn;
    end

    % Collect bond distances if available
    if isfield(atom, 'bond') && ~isempty(atom(i).bond)
        for b = 1:size(atom(i).bond.dist, 1)
            % Bond format is expected to be [neighbor_idx, distance]
            neighbor_idx = atom(i).bond.index(b,2);
            distance = atom(i).bond.dist(b);

            % Add to bond_dists
            if isKey(bond_dists, atom_type)
                bond_dists(atom_type) = [bond_dists(atom_type), distance];
            else
                bond_dists(atom_type) = distance;
            end

            % Collect bond type pairs
            if neighbor_idx <= length(atom) && neighbor_idx > 0 % && isfield(atom, 'type')
                neighbor_type = string(atom(neighbor_idx).type);

                % Use alphabetical ordering of atom types for consistency
                if strcmp(atom_type, neighbor_type) || strcmp(atom_type, '') || strcmp(neighbor_type, '')
                    bond_pair = atom_type;
                elseif strcmp(atom_type, neighbor_type)
                    bond_pair = atom_type;
                else
                    sorted_types = sort([atom_type, neighbor_type]);
                    bond_pair = [sorted_types{1}, '-', sorted_types{2}];
                end

                if isKey(bond_type_pairs, bond_pair)
                    bond_type_pairs(bond_pair) = [bond_type_pairs(bond_pair), distance];
                else
                    bond_type_pairs(bond_pair) = distance;
                end
            end
        end
    end

    % Collect angles if available
    if isfield(atom, 'angle') && ~isempty(atom(i).angle)
        for a = 1:size(atom(i).angle.angle, 1)
            % angle_data format is expected to be [neigh1_idx, neigh2_idx, angle_value]
            neigh1_idx = atom(i).angle.index(a, 1);
            neigh2_idx = atom(i).angle.index(a, 3);
            angle_value = atom(i).angle.angle(a,1);

            % Add to angles
            if isKey(angles, atom_type)
                angles(atom_type) = [angles(atom_type), angle_value];
            else
                angles(atom_type) = angle_value;
            end

            % Special handling for hydrogen-involved angles
            if neigh1_idx <= size(atom,2) && neigh1_idx > 0 % && isfield(atom, 'type')
                neigh1_type = string(atom(neigh1_idx).type);
                if ~isempty(neigh1_type) && length(neigh1_type) >= 1 && strcmpi(neigh1_type(1), 'H')
                    if isKey(h_involved_angles, neigh1_type)
                        h_involved_angles(neigh1_type) = [h_involved_angles(neigh1_type), angle_value];
                    else
                        h_involved_angles(neigh1_type) = angle_value;
                    end
                end
            end

            if neigh2_idx <= length(atom) && neigh2_idx > 0 % && isfield(atom, 'type')
                neigh2_type = string(atom(neigh2_idx).type);
                if ~isempty(neigh2_type) && length(neigh2_type) >= 1 && strcmpi(neigh2_type(1), 'H')
                    if isKey(h_involved_angles, neigh2_type)
                        h_involved_angles(neigh2_type) = [h_involved_angles(neigh2_type), angle_value];
                    else
                        h_involved_angles(neigh2_type) = angle_value;
                    end
                end
            end

            % Collect angle type triplets
            if neigh1_idx <= size(atom,2) && neigh1_idx > 0 && neigh2_idx <= size(atom,2) && neigh2_idx > 0 % && isfield(atom, 'type')
                neigh1_type = string(atom(neigh1_idx).type);
                neigh2_type = string(atom(neigh2_idx).type);

                % Create a unique representation for the angle triplet
                % For A-B-C, where B is the center atom, order A and C alphabetically
                % Handle special cases first (empty strings)
                neigh_type_sorted = sort([neigh1_type,neigh2_type]);

                neigh1_type = neigh_type_sorted(1);
                neigh2_type = neigh_type_sorted(2);

                if isempty(neigh1_type) || isempty(neigh2_type)
                    angle_triplet = strcat(neigh1_type, '-', atom_type, '-', neigh2_type);
                else
                    % Sort the terminal atom types alphabetically
                    if strcmp(neigh1_type, neigh2_type) > 0
                        % neigh2_type comes before neigh1_type alphabetically
                        angle_triplet = strcat(neigh2_type, '-', atom_type, '-', neigh1_type);
                    else
                        % neigh1_type comes before or is equal to neigh2_type
                        angle_triplet = strcat(neigh1_type, '-', atom_type, '-', neigh2_type);
                    end
                end

                % % Create a unique representation for the angle triplet
                % % For A-B-C, where B is the center atom, order as (min(A,C), B, max(A,C))
                % if strcmp(neigh1_type, neigh2_type) || strcmp(neigh1_type, '') || strcmp(neigh2_type, '')
                %     angle_triplet = strcat(neigh1_type, '-', atom_type, '-', neigh2_type);
                % else
                %     if strcmp(neigh1_type, neigh2_type)
                %         angle_triplet = strcat(neigh1_type, '-', atom_type, '-', neigh2_type);
                %     elseif strcmp(neigh1_type, '') || strcmp(neigh2_type, '')
                %         angle_triplet = strcat(neigh1_type, '-', atom_type, '-', neigh2_type);
                %     else
                %         if strcmp(neigh1_type, neigh2_type) < 0
                %             angle_triplet = strcat(neigh1_type, '-', atom_type, '-', neigh2_type);
                %         else
                %             angle_triplet = strcat(neigh2_type, '-', atom_type, '-', neigh1_type);
                %         end
                %     end
                % end


                if isKey(angle_type_triplets, angle_triplet)
                    angle_type_triplets(angle_triplet) = [angle_type_triplets(angle_triplet), angle_value];
                else
                    angle_type_triplets(angle_triplet) = angle_value;
                end
            end
        end
    end

    % Create a sorted neighbor string
    sorted_neighbor_types = sort(neighbor_types);
    all_neighbor_str = '';
    for j = 1:length(sorted_neighbor_types)
        all_neighbor_str = [all_neighbor_str, sorted_neighbor_types{j}];
    end

    charge = 0;
    if isfield(atom, 'charge')
        charge = atom(i).charge;
    end

    all_neighbor_str_list{end+1} = {atom_type, all_neighbor_str, charge};
end

% Create a consolidated dictionary of unique atom types, their coordination, and charges
unique_patterns = containers.Map; %('KeyType', 'char', 'ValueType', 'any');
for i = 1:length(all_neighbor_str_list)
    atom_info = all_neighbor_str_list{i};
    atom_type = atom_info{1};
    neighbor_str = atom_info{2};
    charge = atom_info{3};

    key = [atom_type, '_', neighbor_str];
    key = strcat(atom_type, '_', neighbor_str);

    if isKey(unique_patterns, key)
        data = unique_patterns(key);
        data.count = data.count + 1;
        data.charges = [data.charges, charge];
        unique_patterns(key) = data;
    else
        data = struct();
        data.count = 1;
        data.charges = charge;
        data.atom_type = atom_type;
        data.neighbor_pattern = neighbor_str;
        unique_patterns(key) = data;
    end
end

% Format output
output = {};

% Calculate total mass for density calculation
total_mass = 0;
if isfield(atom, 'mass')
    for i = 1:length(atom)
        total_mass = total_mass + atom(i).mass;
    end
end

% Add box dimensions, volume, and density information
if ~isempty(Box_dim) || ~isempty(Cell)
    output{end+1} = 'System Dimensions and Properties';
    output{end+1} = repmat('-', 1, 80);

    % Box_dim representation
    if ~isempty(Box_dim)
        output{end+1} = 'Box_dim (Å):';
        if length(Box_dim) == 3
            output{end+1} = sprintf('  Orthogonal box: [%.4f, %.4f, %.4f]', Box_dim(1), Box_dim(2), Box_dim(3));
            volume = Box_dim(1) * Box_dim(2) * Box_dim(3);
            output{end+1} = sprintf('  Volume: %.4f Å³', volume);
        elseif length(Box_dim) == 6
            output{end+1} = sprintf('  Triclinic box: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]', ...
                Box_dim(1), Box_dim(2), Box_dim(3), Box_dim(4), Box_dim(5), Box_dim(6));
        elseif length(Box_dim) == 9
            output{end+1} = sprintf('  Full matrix: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]', ...
                Box_dim(1), Box_dim(2), Box_dim(3), Box_dim(4), Box_dim(5), Box_dim(6), ...
                Box_dim(7), Box_dim(8), Box_dim(9));
        end
    end

    % Cell representation
    if ~isempty(Cell)
        a = Cell(1);
        b = Cell(2);
        c = Cell(3);
        alpha = Cell(4);
        beta = Cell(5);
        gamma = Cell(6);

        output{end+1} = 'Cell parameters:';
        output{end+1} = sprintf('  a, b, c (Å): %.4f, %.4f, %.4f', a, b, c);
        output{end+1} = sprintf('  α, β, γ (°): %.4f, %.4f, %.4f', alpha, beta, gamma);

        % Calculate volume from cell parameters
        if abs(alpha - 90) < 1e-6 && abs(beta - 90) < 1e-6 && abs(gamma - 90) < 1e-6
            % Orthogonal box
            volume = a * b * c;
        else
            % Triclinic box, use the general formula
            alpha_rad = alpha * pi / 180;
            beta_rad = beta * pi / 180;
            gamma_rad = gamma * pi / 180;
            volume = a * b * c * sqrt(1 - cos(alpha_rad)^2 - cos(beta_rad)^2 - ...
                cos(gamma_rad)^2 + 2 * cos(alpha_rad) * cos(beta_rad) * cos(gamma_rad));
        end
        output{end+1} = sprintf('  Volume: %.4f Å³', volume);
    end

    % Calculate and display density
    % Only calculate if we have a volume and total mass is non-zero
    if exist('volume', 'var') && total_mass > 0
        % Convert from amu/Å³ to g/cm³
        density = total_mass / volume * AMU_TO_G_PER_CM3;
        output{end+1} = 'System properties:';
        output{end+1} = sprintf('  Total mass: %.4f amu', total_mass);
        output{end+1} = sprintf('  Density: %.4f g/cm³', density);
    end

    % Explanation of variables
    output{end+1} = '';
    output{end+1} = 'Box_dim and Cell explanations:';
    output{end+1} = '  Box_dim: A 1D array of box dimensions, typically in Angstroms.';
    output{end+1} = '          For orthogonal boxes: [Lx, Ly, Lz]';
    output{end+1} = '          For triclinic boxes: [Lx, Ly, Lz, xy, xz, yz] or [Lx, Ly, Lz, α, β, γ]';
    output{end+1} = '  Cell: A 1×6 array with cell parameters [a, b, c, α, β, γ]';
    output{end+1} = '        where a, b, c are lengths and α, β, γ are angles in degrees.';
    output{end+1} = repmat('-', 1, 80);
    output{end+1} = '';
end

if isfield(atom,'charge')
    total_charge=sum([atom.charge]);
else
    total_charge=0;
end
output{end+1} = sprintf('Total charge: %.7f', total_charge);

if abs(round(total_charge) - total_charge) > 1e-10
    output{end+1} = 'Warning: Non-integer total charge. Adjusting to nearest integer.';
    target_charge = round(total_charge);

    % Calculate sum of charges
    sum_charges = 0;
    if isfield(atom, 'charge')
        for i = 1:length(atom)
            sum_charges = sum_charges + atom(i).charge;
        end
    end

    output{end+1} = sprintf('Final total charge: %.7f (target was %d)', sum_charges, target_charge);
end

output{end+1} = '';
output{end+1} = 'Unique Atom Types and Their Coordination Environment';
output{end+1} = repmat('-', 1, 80);
output{end+1} = sprintf('%-10s %-6s %-20s %15s', 'Type', 'Count', 'Neighbors', 'Charge');
output{end+1} = repmat('-', 1, 80);

% Sort by atom type for a more organized display
keys = unique_patterns.keys;
sorted_keys = sort(keys);

for i = 1:length(sorted_keys)
    key = sorted_keys{i};
    data = unique_patterns(key);
    atom_type = data.atom_type;
    neighbor_pattern = data.neighbor_pattern;
    count = data.count;
    charges = data.charges;

    % Find unique charge values (with some tolerance for floating point comparison)
    if ~isempty(charges)
        if length(charges) == 1
            unique_charges = charges;
        else
            unique_charges = [];
            for j = 1:length(charges)
                charge = charges(j);
                is_unique = true;

                for k = 1:length(unique_charges)
                    if abs(charge - unique_charges(k)) < 1e-6
                        is_unique = false;
                        break;
                    end
                end

                if is_unique
                    unique_charges = [unique_charges, charge];
                end
            end
        end

        % If there's only one unique charge, display it as before
        if length(unique_charges) == 1
            charge_str = sprintf('%.5f', unique_charges(1));
        else
            % Otherwise, display all unique charges separated by commas
            charge_str = '';
            sorted_unique_charges = sort(unique_charges);
            for j = 1:length(sorted_unique_charges)
                if j > 1
                    charge_str = [charge_str, ', '];
                end
                charge_str = [charge_str, sprintf('%.5f', sorted_unique_charges(j))];
            end
        end

        output{end+1} = sprintf('%-10s %-6d %-20s %15s', atom_type, count, neighbor_pattern, charge_str);
    else
        output{end+1} = sprintf('%-10s %-6d %-20s %15s', atom_type, count, neighbor_pattern, 'N/A');
    end
end

output{end+1} = repmat('-', 1, 80);

% Add detailed statistics for average coordination, bond distances, and angles
output{end+1} = '';
output{end+1} = 'Detailed Atom Type Statistics with Standard Deviations';
output{end+1} = repmat('-', 1, 80);
header = sprintf('%-10s %-6s %-16s %-20s %-20s', 'Type', 'Count', 'Coord#', 'Bond Dist (Å)', 'Angle (°)');
output{end+1} = header;
output{end+1} = repmat('-', 1, 80);

atom_types = atom_type_counts.keys;
sorted_atom_types = sort(atom_types);

for i = 1:length(sorted_atom_types)
    atom_type = sorted_atom_types{i};
    count = atom_type_counts(atom_type);

    % Calculate average coordination number and std dev
    if isKey(coord_nums, atom_type)
        cn_data = coord_nums(atom_type);
        if length(cn_data) > 0
            avg_cn = mean(cn_data);
            if length(cn_data) > 1
                std_cn = std(cn_data);
            else
                std_cn = 0;
            end
            cn_str = sprintf('%.2f ± %.2f', avg_cn, std_cn);
        else
            cn_str = 'N/A';
        end
    else
        cn_str = 'N/A';
    end

    % Calculate average bond distance and std dev
    if isKey(bond_dists, atom_type)
        dist_data = bond_dists(atom_type);
        if length(dist_data) > 0
            avg_dist = mean(dist_data);
            if length(dist_data) > 1
                std_dist = std(dist_data);
            else
                std_dist = 0;
            end
            dist_str = sprintf('%.4f ± %.4f', avg_dist, std_dist);
        else
            dist_str = 'N/A';
        end
    else
        dist_str = 'N/A';
    end

    % Calculate average angle and std dev
    angle_str = 'N/A';
    if isKey(angles, atom_type)
        angle_data = angles(atom_type);
        if length(angle_data) > 0
            avg_angle = mean(angle_data);
            if length(angle_data) > 1
                std_angle = std(angle_data);
            else
                std_angle = 0;
            end
            angle_str = sprintf('%.3f ± %.3f', avg_angle, std_angle);
        end
    elseif ~isempty(atom_type) && length(atom_type) >= 1 && strcmpi(atom_type(1), 'H') && isKey(h_involved_angles, atom_type)
        % For hydrogen atom types with no angles (as center atom), use the angles where H is involved
        h_angle_data = h_involved_angles(atom_type);
        if length(h_angle_data) > 0
            avg_angle = mean(h_angle_data);
            if length(h_angle_data) > 1
                std_angle = std(h_angle_data);
            else
                std_angle = 0;
            end
            angle_str = sprintf('%.3f ± %.3f (H-O-M)', avg_angle, std_angle);
        end
    end

    % Format the output line
    output{end+1} = sprintf('%-10s %-6d %-16s %-20s %-20s', atom_type, count, cn_str, dist_str, angle_str);
end

output{end+1} = repmat('-', 1, 80);

% Add bond statistics for unique atom type pairs
output{end+1} = '';
output{end+1} = 'Bond Statistics for Unique Atom Type Pairs';
output{end+1} = repmat('-', 1, 80);
output{end+1} = sprintf('%-25s %-8s %-20s', 'Bond Pair', 'Count', 'Distance (Å)');
output{end+1} = repmat('-', 1, 80);

% Sort by bond pair for a more organized display
bond_pairs = bond_type_pairs.keys;
sorted_bond_pairs = sort(bond_pairs);

for i = 1:length(sorted_bond_pairs)
    bond_pair = sorted_bond_pairs{i};
    distances = bond_type_pairs(bond_pair);
    count = length(distances);

    if count > 0
        avg_dist = mean(distances);
        if count > 1
            std_dist = std(distances);
        else
            std_dist = 0;
        end

        output{end+1} = sprintf('%-25s %-8d %.4f ± %.4f', bond_pair, count, avg_dist, std_dist);
    end
end

output{end+1} = repmat('-', 1, 80);

% Add angle statistics for unique atom type triplets
output{end+1} = '';
output{end+1} = 'Angle Statistics for Unique Atom Type Triplets';
output{end+1} = repmat('-', 1, 80);
output{end+1} = sprintf('%-25s %-8s %-20s', 'Angle Triplet', 'Count', 'Angle (°)');
output{end+1} = repmat('-', 1, 80);

% Sort by angle triplet for a more organized display
angle_triplets = angle_type_triplets.keys;
sorted_angle_triplets = sort(angle_triplets);

for i = 1:length(sorted_angle_triplets)
    angle_triplet = sorted_angle_triplets{i};
    angles_data = angle_type_triplets(angle_triplet);
    count = length(angles_data);

    if count > 0
        avg_angle = mean(angles_data);
        if count > 1
            std_angle = std(angles_data);
        else
            std_angle = 0;
        end

        output{end+1} = sprintf('%-25s %-8d %.3f ± %.3f', angle_triplet, count, avg_angle, std_angle);
    end
end

output{end+1} = repmat('-', 1, 80);

% Join everything into a string
result = sprintf('%s\n', output{:});

% Write to log file if specified
if ~isempty(log_file)
    fid = fopen(log_file, 'w');
    if fid ~= -1
        fprintf(fid, '%s\n', result);
        fclose(fid);
    else
        warning('Could not open log file %s for writing', log_file);
    end
end

end
