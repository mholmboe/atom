%% cell_list_dist_matrix_atom.m
% * This function calculates the distance matrix from the atom struct for
% a triclinic box and PBC using a cell-list approach. Additionally, it 
% returns a full NxN distance matrix and three NxN displacement matrices 
% for pairs within the cutoff.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # [dist_matrix,bond_list,dist_list,X_dist,Y_dist,Z_dist] = cell_list_dist_matrix_atom(atom,Box_dim) % Basic input arguments
% # [dist_matrix,bond_list,dist_list,X_dist,Y_dist,Z_dist] = cell_list_dist_matrix_atom(atom,Box_dim,1.25) % Setting the max distance for bonds with H's
% # [dist_matrix,bond_list,dist_list,X_dist,Y_dist,Z_dist] = cell_list_dist_matrix_atom(atom,Box_dim,1.25,2.25) % Setting the max distance for bonds, otherwise default is 12 Å
% # [dist_matrix,bond_list,dist_list,X_dist,Y_dist,Z_dist] = cell_list_dist_matrix_atom(atom,Box_dim,1.25,2.25,2.25) % Setting the max distance for angles, default is not to calculate any angles

function [dist_matrix,bond_list,dist_list,X_dist,Y_dist,Z_dist] = cell_list_dist_matrix_atom(atom, Box_dim,varargin)

% Parse input & setup
if nargin==2
    rmaxlong = 2.45;
elseif nargin==3
    rmaxlong = varargin{1};
elseif nargin==4
    rmaxshort = varargin{1};
    rmaxlong = varargin{2};
end

if nargin==5
    rAtomType=varargin{3};
else
    rAtomType='H';
end

% Possibly convert Box_dim into [a,b,c,alpha,beta,gamma] form
if numel(Box_dim)==9
    Cell = Box_dim2Cell(Box_dim);
elseif numel(Box_dim)==6
    Cell = Box_dim;
elseif numel(Box_dim)==3
    Cell = Box_dim;   % Orthogonal box with angles=90
end

a = Cell(1);
b = Cell(2);
c = Cell(3);

% For orthogonal, angles = 90
if numel(Cell)==3
    alpha_rad = deg2rad(90);
    beta_rad  = deg2rad(90);
    gamma_rad = deg2rad(90);
else
    alpha_rad = deg2rad(Cell(4));
    beta_rad  = deg2rad(Cell(5));
    gamma_rad = deg2rad(Cell(6));
end

% Number of atoms
N = numel(atom);

% Extract coordinates into Nx3 matrix
coords = [[atom.x]' [atom.y]' [atom.z]'];

% Construct triclinic box matrix H
ax = a;
bx = b*cos(gamma_rad);
by = b*sin(gamma_rad);
cx = c*cos(beta_rad);
cy = c*(cos(alpha_rad) - cos(beta_rad)*cos(gamma_rad)) / sin(gamma_rad);
cz = sqrt( c^2 - cx^2 - cy^2 );

H = [ax bx cx;
    0  by cy;
    0   0  cz];
Hinv = inv(H);

% Build the cell list
%   (Pick a cell size that ensures neighbors can be found)
if exist('rmaxlong','var')
    cellSize = 2*rmaxlong;
else
    cellSize = 5;  % fallback
end

boxCart = H * [1; 1; 1]; % bounding box corner
boundingBoxSize = [max(abs([ax bx cx])), ...
    max(abs([0  by cy])), ...
    max(abs([0   0  cz]))];
nCells = max(floor(boundingBoxSize / cellSize), 1);

% Convert coords -> fractional
fracCoords = (Hinv * coords')';
fracCoords = fracCoords - floor(fracCoords);  % keep in [0,1)

% Compute which cell each atom belongs to
cellIndex = floor(fracCoords .* nCells);
% fix boundaries
outOfBound = (cellIndex(:,1) >= nCells(1));
cellIndex(outOfBound,1) = nCells(1) - 1;
outOfBound = (cellIndex(:,2) >= nCells(2));
cellIndex(outOfBound,2) = nCells(2) - 1;
outOfBound = (cellIndex(:,3) >= nCells(3));
cellIndex(outOfBound,3) = nCells(3) - 1;

% Convert (ix,iy,iz) -> single linear index
cellLinIdx = sub2ind(nCells, ...
    cellIndex(:,1)+1, ...
    cellIndex(:,2)+1, ...
    cellIndex(:,3)+1);
numCells = prod(nCells);
cellList = cell(numCells,1);

for iAtom = 1:N
    cIdx = cellLinIdx(iAtom);
    cellList{cIdx} = [cellList{cIdx}, iAtom];
end

% Identify neighbor cells (27 total with PBC)
neighborOffsets = -1:1;
neighborIndices = zeros(27,3);
count = 0;
for ix = neighborOffsets
    for iy = neighborOffsets
        for iz = neighborOffsets
            count = count + 1;
            neighborIndices(count,:) = [ix, iy, iz];
        end
    end
end

% Loop over cells/neighbors & compute distances
bond_list = [];
dist_list = [];

% Initialize NxN output matrices
% (Set them to 0 for now; distances for i=j or outside cutoff remain 0)
dist_matrix = zeros(N, N);
X_dist   = zeros(N, N);
Y_dist   = zeros(N, N);
Z_dist   = zeros(N, N);

for cID = 1:numCells
    atomListC = cellList{cID};
    if isempty(atomListC), continue; end

    [cx, cy, cz] = ind2sub(nCells, cID);

    for nIdx = 1:size(neighborIndices,1)
        dxIdx = neighborIndices(nIdx,1);
        dyIdx = neighborIndices(nIdx,2);
        dzIdx = neighborIndices(nIdx,3);

        nx = mod(cx-1+dxIdx, nCells(1)) + 1;
        ny = mod(cy-1+dyIdx, nCells(2)) + 1;
        nz = mod(cz-1+dzIdx, nCells(3)) + 1;

        neighborCellLin = sub2ind(nCells, nx, ny, nz);
        atomListN = cellList{neighborCellLin};
        if isempty(atomListN), continue; end

        % Avoid double-counting
        if (neighborCellLin < cID)
            continue;
        end

        %%
        if exist('rmaxshort','var')
            % Compute pairwise distances iAtom <-> jAtom
            for i = 1:length(atomListC)
                iAtom = atomListC(i);
                for j = 1:length(atomListN)
                    jAtom = atomListN(j);

                    % If in the same cell, only do jAtom>iAtom
                    if (neighborCellLin == cID && jAtom <= iAtom)
                        continue;
                    end

                    % Check hydrogen vs. non-hydrogen
                    isH_i = strcmp(atom(iAtom).type,rAtomType);
                    isH_j = strcmp(atom(jAtom).type,rAtomType);

                    if isH_i || isH_j
                        localCutoff = rmaxshort; % short cutoff for hydrogen
                    else
                        localCutoff = rmaxlong;  % default cutoff
                    end

                    % Minimum-image in fractional space
                    diffFrac = (Hinv*(coords(jAtom,:)' - coords(iAtom,:)'))';
                    diffFrac = diffFrac - round2dec(diffFrac);  % wrap to [-0.5,0.5)

                    % Convert back to Cartesian
                    dVec = (H * diffFrac')';
                    d = norm(dVec);

                    if d <= localCutoff %cutoff
                        % store pair & distance
                        bond_list = [bond_list; iAtom, jAtom];
                        dist_list = [dist_list; d];

                        % fill NxN distance matrix
                        dist_matrix(iAtom, jAtom) = d;
                        dist_matrix(jAtom, iAtom) = d;  % symmetric

                        % fill NxN displacement matrices
                        X_dist(iAtom, jAtom) = dVec(1);
                        X_dist(jAtom, iAtom) = -dVec(1);

                        Y_dist(iAtom, jAtom) = dVec(2);
                        Y_dist(jAtom, iAtom) = -dVec(2);

                        Z_dist(iAtom, jAtom) = dVec(3);
                        Z_dist(jAtom, iAtom) = -dVec(3);
                    end
                end
            end
            %%
        else
            % Compute pairwise distances iAtom <-> jAtom
            for i = 1:length(atomListC)
                iAtom = atomListC(i);
                for j = 1:length(atomListN)
                    jAtom = atomListN(j);

                    % If in the same cell, only do jAtom>iAtom
                    if (neighborCellLin == cID && jAtom <= iAtom)
                        continue;
                    end

                    % Minimum-image in fractional space
                    diffFrac = (Hinv*(coords(jAtom,:)' - coords(iAtom,:)'))';
                    diffFrac = diffFrac - round2dec(diffFrac);  % wrap to [-0.5,0.5)

                    % Convert back to Cartesian
                    dVec = (H * diffFrac')';
                    d = norm(dVec);

                    if d <= rmaxlong
                        % store pair & distance
                        bond_list = [bond_list; iAtom, jAtom];
                        dist_list = [dist_list; d];

                        % fill NxN distance matrix
                        dist_matrix(iAtom, jAtom) = d;
                        dist_matrix(jAtom, iAtom) = d;  % symmetric

                        % fill NxN displacement matrices
                        X_dist(iAtom, jAtom) = dVec(1);
                        X_dist(jAtom, iAtom) = -dVec(1);

                        Y_dist(iAtom, jAtom) = dVec(2);
                        Y_dist(jAtom, iAtom) = -dVec(2);

                        Z_dist(iAtom, jAtom) = dVec(3);
                        Z_dist(jAtom, iAtom) = -dVec(3);
                    end
                end
            end
        end
    end
end

assignin('caller','pair',bond_list);
assignin('caller','dist',dist_list);
assignin('caller','X_dist',(X_dist)');
assignin('caller','Y_dist',(Y_dist)');
assignin('caller','Z_dist',(Z_dist)');

end
