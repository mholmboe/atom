%% import_atom_gro.m
% * This function import .gro files into an atom struct variable
% * varargin can be used to translate, alt. center+translate the molecule
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom_gro('molecule.gro')
% # atom = import_atom_gro('molecule.gro',[10 5 2])
% # atom = import_atom_gro('molecule.gro',[10 5 0],[35.24 24.23 52.23])


function [atom, Box_dim] = import_atom_gro(filename, varargin)

fileID = fopen(filename, 'r');

% Read header lines (Title and nAtoms)
fgetl(fileID); % Skip Title line
Line2 = fgetl(fileID); % Read the number of atoms
nAtoms = str2double(Line2);

% Preallocate arrays for atom attributes
MolID = zeros(nAtoms, 1);
Resname = cell(nAtoms, 1);
AtomType = cell(nAtoms, 1);
X_coord = zeros(nAtoms, 1);
Y_coord = zeros(nAtoms, 1);
Z_coord = zeros(nAtoms, 1);
VeloX = NaN(nAtoms, 1);  % Default NaN for velocities if missing
VeloY = NaN(nAtoms, 1);
VeloZ = NaN(nAtoms, 1);

% Initialize atom structure fields
atom = struct('molid', {}, 'resname', {}, 'type', {}, 'fftype', {}, 'index', [], ...
    'neigh', struct('type', {{}}, 'index', zeros(6, 1), 'dist', zeros(6, 1)), ...
    'bond', struct('type', zeros(6, 1), 'index', zeros(6, 1)), ...
    'angle', struct('type', zeros(6, 1), 'index', zeros(6, 1)), ...
    'x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', []);

nmol = 1; first_in = zeros(nAtoms, 1); last_in = zeros(nAtoms, 1);
for i = 1:nAtoms
    line = fgetl(fileID);

    % Parse the current line using string indexing and convert to appropriate types
    MolID(i) = str2double(line(1:5));
    Resname{i} = strtrim(line(6:10));
    AtomType{i} = strtrim(line(11:15));
    X_coord(i) = 10 * str2double(line(21:28));
    Y_coord(i) = 10 * str2double(line(29:36));
    Z_coord(i) = 10 * str2double(line(37:44));

    % Check if velocity is available, otherwise assign NaN
    if numel(line) > 45
        VeloX(i) = str2double(line(45:52));
        VeloY(i) = str2double(line(53:60));
        VeloZ(i) = str2double(line(61:68));
    end

    % Assign molecule ID and structure
    if i > 1 && MolID(i) ~= MolID(i - 1)
        nmol = nmol + 1;
        atom(i).molid = nmol;
        first_in(nmol) = i;
        last_in(nmol - 1) = i - 1;
    else
        atom(i).molid = nmol;
    end

    % Update atom struct with parsed data
    atom(i).resname = Resname(i);
    atom(i).type = AtomType(i);
    atom(i).fftype = AtomType(i);  % Assuming fftype is the same as type
    atom(i).index = mod(i, 100000);
    atom(i).neigh.index = [0;0;0;0;0;0];
    atom(i).neigh.dist = [0;0;0;0;0;0];
    atom(i).bond.type = [0;0;0;0;0;0];
    atom(i).bond.index = [0;0;0;0;0;0];
    atom(i).angle.type = [0;0;0;0;0;0];
    atom(i).angle.index = [0;0;0;0;0;0];
    atom(i).x = X_coord(i);
    atom(i).y = Y_coord(i);
    atom(i).z = Z_coord(i);
    atom(i).vx = VeloX(i);
    atom(i).vy = VeloY(i);
    atom(i).vz = VeloZ(i);
end
last_in(nmol) = nAtoms;  % Handle the last molecule

Box_string = fgetl(fileID);
fclose(fileID);

Box_dim = str2double(strsplit(char(Box_string))) * 10;
Box_dim(isnan(Box_dim)) = [];  % Clean up NaN values

% If optional translation or centering is specified, apply it
if nargin >= 2
    atom = translate_atom(atom, varargin{1} + [0 0 -median([atom.z])], 'all');
end

if nargin >= 3
    atom = center_atom(atom, varargin{2}, 'all', 'xyz');
    atom = translate_atom(atom, varargin{1} + [0 0 -median([atom.z])], 'all');
end

% Prepare XYZ data and labels for output
XYZ_data = [[atom.x]' [atom.y]' [atom.z]'];
XYZ_labels = {atom.type}';

Cell = Box_dim2Cell(Box_dim);

% Assign variables to the caller workspace
assignin('caller', 'XYZ_labels', XYZ_labels);
assignin('caller', 'XYZ_data', XYZ_data);
assignin('caller', 'atom', atom);
assignin('caller', 'nAtoms', nAtoms);
assignin('caller', 'Box_dim', Box_dim);
assignin('caller', 'Cell', Cell);
assignin('caller','MolID',MolID)

disp('.gro file imported');

end
