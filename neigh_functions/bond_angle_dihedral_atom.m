%% bond_angle_dihedral_atom.m
% * This function tries to find all bonds, angles and the dihedrals of the atom struct.
% * Compute bonds, angles, dihedrals and impropers based on atom struct.
% * Uses bond_atom() to generate Bond_index and Angle_index, then
% * derives Dihedral_index from the bond graph without duplicates.
% * Box_dim is the box dimension vector
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=bond_angle_dihedral_atom(atom,Box_dim) % Basic input arguments
% # atom=bond_angle_dihedral_atom(atom) % When the PBC is not important
% # atom=bond_angle_dihedral_atom(atom,Box_dim,1.25,2.25) % Setting the max distance rmaxshort and rmaxlong for bonds with H's
% # atom=bond_angle_dihedral_atom(atom,Box_dim,1.25,2.25,'more') % Will write more info to the calling workspace

function atom = bond_angle_dihedral_atom(atom, varargin)


    if size(atom,2) > 10000
        disp('This is a large molecule or system, are you sure you want to calculate all dihedrals?')
        disp('If not, use the bond_atom() or the bond_angle_atom() functions!')
        pause(2)
    end

    disp('Calculating bonds and angles')

    % -----------------------------
    %  Handle input arguments
    % -----------------------------
    if nargin < 2
        % Dummy Box_dim when PBC is not important
        Box_dim = 1e6 * [1 1 1];
        rmaxshort = 1.25;
        rmaxlong  = 2.25;
    elseif nargin == 2
        Box_dim   = varargin{1};
        rmaxshort = 1.25;
        rmaxlong  = 2.25;
    elseif nargin == 3
        Box_dim   = varargin{1};
        rmaxshort = varargin{2};
        rmaxlong  = 2.25;
    else
        Box_dim   = varargin{1};
        rmaxshort = varargin{2};
        rmaxlong  = varargin{3};
    end

    % ---------------------------------------------
    % Call bond_atom to generate Bond_index, etc.
    % (Assumes bond_atom assigns Bond_index, Angle_index, dist_matrix, etc.
    %  in this function's workspace via assignin('caller',...))
    % ---------------------------------------------
    atom = bond_atom(atom, Box_dim, rmaxlong);

    % At this point, Bond_index and Angle_index should exist in this workspace
    if ~exist('Bond_index','var') || isempty(Bond_index)
        warning('Bond_index not found or empty after bond_atom. No dihedrals will be generated.');
        Dihedral_index          = [];
        Improper_dihedral_index = [];
        Pairlist                = [];
        nBonds                  = 0;
        nAngles                 = 0;
        nDihedrals              = 0;
        nImpropers              = 0;
        dist_matrix             = [];
        % Push variables to caller and return
        assignin('caller','Bond_index',Bond_index);
        assignin('caller','Angle_index',Angle_index);
        assignin('caller','Dihedral_index',Dihedral_index);
        assignin('caller','Improper_dihedral_index',Improper_dihedral_index);
        assignin('caller','Pairlist',Pairlist);
        assignin('caller','nBonds',nBonds);
        assignin('caller','nAngles',nAngles);
        assignin('caller','nDihedrals',nDihedrals);
        assignin('caller','nImpropers',nImpropers);
        assignin('caller','dist_matrix',dist_matrix);
        return
    end

    % If Angle_index not produced by bond_atom, make sure it exists
    if ~exist('Angle_index','var')
        Angle_index = [];
    end

    % ---------------------------------------------
    % Basic counters (if bond_atom set them; otherwise infer)
    % ---------------------------------------------
    if exist('nBonds','var')
        nBonds = nBonds;
    else
        nBonds = size(Bond_index,1);
    end

    if exist('nAngles','var')
        nAngles = nAngles;
    else
        nAngles = size(Angle_index,1);
    end

    % ---------------------------------------------
    % Build dihedrals from bond graph (no duplicates)
    % ---------------------------------------------
    disp('Calculating dihedrals and pairlist')

    % Ensure bond list is sorted with smaller index first (for consistency)
    Bond_index(:,1:2) = sort(Bond_index(:,1:2), 2);

    % Total number of atoms (from bonds)
    N_atom = max(Bond_index(:));

    % Build adjacency list
    adjacency_list = cell(N_atom, 1);
    for i = 1:size(Bond_index, 1)
        a1 = Bond_index(i, 1);
        a2 = Bond_index(i, 2);
        adjacency_list{a1} = [adjacency_list{a1}, a2];
        adjacency_list{a2} = [adjacency_list{a2}, a1];
    end

    % Loop over bonds and make dihedrals A-B-C-D where B-C is a bond
    Dihedral_index = [];
    for i = 1:size(Bond_index, 1)
        atom2 = Bond_index(i, 1);  % B
        atom3 = Bond_index(i, 2);  % C

        % Neighbors of atom2 excluding atom3
        neighbors2 = adjacency_list{atom2};
        neighbors2(neighbors2 == atom3) = [];

        % Neighbors of atom3 excluding atom2
        neighbors3 = adjacency_list{atom3};
        neighbors3(neighbors3 == atom2) = [];

        for atom1 = neighbors2     % A
            for atom4 = neighbors3 % D
                if atom1 == atom4
                    % Skip 3-member cycles A-B-C-A
                    continue
                end

                dih = [atom1, atom2, atom3, atom4];

                % Enforce canonical ordering to avoid A-B-C-D vs D-C-B-A
                if atom1 < atom4
                    Dihedral_index = [Dihedral_index; dih];
                else
                    Dihedral_index = [Dihedral_index; dih([4 3 2 1])];
                end
            end
        end
    end

    % Remove duplicate dihedrals and sort
    if ~isempty(Dihedral_index)
        Dihedral_index = unique(Dihedral_index, 'rows');
        % Sort primarily by atom2, secondarily by atom3 (column-wise priority)
        Dihedral_index = sortrows(Dihedral_index, [2 3]);
        nDihedrals = size(Dihedral_index, 1);
    else
        nDihedrals = 0;
    end

    % Pairlist from dihedral terminal atoms, excluding existing bonds
    if ~isempty(Dihedral_index)
        Pairlist = Dihedral_index(:,[1 4]);
        Pairlist = unique(Pairlist, 'rows');
        % Remove any pairs that are already bonds
        Pairlist(ismember(Pairlist, Bond_index(:,1:2), 'rows'), :) = [];
    else
        Pairlist = [];
    end

    % ---------------------------------------------
    % Improper dihedrals (planarity-based) – unchanged logic,
    % but kept separate from Dihedral_index
    % ---------------------------------------------
    disp('Calculating improper dihedrals with planarity check...');

    Improper_dihedral_index = [];
    planarityThreshold = 1e-1;

    % Only proceed if we have Angle_index information
    if ~isempty(Angle_index)
        for central_atom = 1:N_atom

            ind3 = find(Angle_index(:,2) == central_atom);

            % Example condition: exactly three angles with this central atom
            % and specific atom type (adapt as needed)
            if numel(ind3) == 3 && ...
               sum(ismember(atom(central_atom).type, {'CG2O1','CG2O3'})) > 0

                % Collect coordinates from the three angles
                pos_vec = [Angle_index(ind3(1),5:7); ...
                           Angle_index(ind3(1),8:10); ...
                           Angle_index(ind3(2),5:7); ...
                           Angle_index(ind3(2),8:10); ...
                           Angle_index(ind3(3),5:7); ...
                           Angle_index(ind3(3),8:10)];
                pos_vec = unique(pos_vec, 'rows');

                if size(pos_vec,1) < 3
                    continue
                end

                vec1 = pos_vec(1,:);
                vec2 = pos_vec(2,:);
                vec3 = pos_vec(3,:);

                % Compute normals to check planarity
                normal1 = cross(vec1, vec2);
                normal2 = cross(vec1, vec3);

                if norm(normal1) < 1e-12 || norm(normal2) < 1e-12
                    continue
                end

                planarity = dot(normal1, normal2) / (norm(normal1)*norm(normal2));

                if abs(abs(planarity) - 1) < planarityThreshold
                    % Derive the three neighbor atom indices from Angle_index
                    neighAtoms = unique(Angle_index(ind3, [1 3]));
                    dihedralRow = [central_atom, neighAtoms(:)'];
                    Improper_dihedral_index = [Improper_dihedral_index; dihedralRow];
                end
            end
        end

        if ~isempty(Improper_dihedral_index)
            Improper_dihedral_index = unique(Improper_dihedral_index, 'rows');
            Improper_dihedral_index = sortrows(Improper_dihedral_index);
        end
    end

    nImpropers = size(Improper_dihedral_index, 1);

    % ---------------------------------------------
    % Export variables to caller (as in your original)
    % ---------------------------------------------
    if ~exist('dist_matrix','var')
        dist_matrix = [];
    end

    assignin('caller','dist_matrix',dist_matrix);
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','Angle_index',Angle_index);
    assignin('caller','nBonds',nBonds);
    assignin('caller','nAngles',nAngles);
    assignin('caller','nDihedrals',nDihedrals);
    assignin('caller','nImpropers',nImpropers);
    assignin('caller','Dihedral_index',Dihedral_index);
    assignin('caller','Improper_dihedral_index',Improper_dihedral_index);
    assignin('caller','Pairlist',Pairlist);
end
