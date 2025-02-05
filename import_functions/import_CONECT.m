%% import_CONECT
% * % Script to extract CONECT records from PDB file and store them in the
% Bond_index variable. If the .pdb file alos contains Box_dim info, the
% bonddistances are also calculated.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # Bond_index = import_CONECT(filename) % Basic input arguments
%
function Bond_index = import_CONECT(filename)

try
    [atom,Box_dim]=import_atom_pdb(filename);
catch
    disp('Could not import the atom struct or Box_dim');
end

% Open file for reading
fid = fopen(filename, 'r');

% Initialize Bond_index
Bond_index = [];

% Check if the file opened successfully
if fid == -1
    error('Could not open file %s', filename);
end

% Read the file line by line
while ~feof(fid)
    tline = fgetl(fid); % Read one line
    if startsWith(tline, 'CONECT') % Process only CONECT records
        % Extract indices
        indices = sscanf(tline(7:end), '%d'); % Read all numbers after 'CONECT'
        if numel(indices) > 1
            atom_index = indices(1); % First index (central atom)
            bonded_indices = indices(2:end); % Bonded atoms
            % Append to Bond_index
            for j = 1:length(bonded_indices)
                Bond_index = [Bond_index; atom_index, bonded_indices(j)];
            end
        end
    end
end

% Close the file
fclose(fid);

% Post-process Bond_index to ensure smallest index in first column
Bond_index = sort(Bond_index, 2); % Sort each row so the smaller index is first

% Remove duplicate pairs
Bond_index = unique(Bond_index, 'rows');

if exist("atom","var")
assignin('caller',"Bond_index_CONECT",Bond_index);
    Bond_index=[Bond_index Bond_index(:,1)];
    atom=recalc_bond_atom(atom,Box_dim,Bond_index);

end
