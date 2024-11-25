%% bond_angle_dihedral_atom.m
% * This function tries to find all bonds, angles and the dihedrals of the atom struct.
% *
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

function atom = bond_angle_dihedral_atom(atom,varargin)

if size(atom,2)>10000
    disp('This is a large molecule or system, are you sure you want to calculate all dihedrals?')
    disp('If not, use the bond_atom() or the bond_angle_atom() functions!')
    pause(2)
end

disp('Calculating bonds and angles')

if nargin<=4
    if nargin<4
        if nargin==1
            Box_dim=1e6*[1 1 1]; % Dummy Box_dim, when the PBC is not important
        else
            Box_dim=varargin{1};
        end
        rmaxshort=1.25;
        rmaxlong=2.25;
    else
        Box_dim=varargin{1};
        rmaxshort=varargin{2};
        rmaxlong=varargin{3};
    end
    atom=bond_angle_atom(atom,Box_dim,rmaxshort,rmaxlong);
elseif nargin>4
    Box_dim=varargin{1};
    rmaxshort=varargin{2};
    rmaxlong=varargin{3};
    atom=bond_angle_atom(atom,Box_dim,rmaxshort,rmaxlong,'more');
end

Dihedral_index=[];Pairlist=[];
nDihedrals=size(Dihedral_index,2);
if size(Bond_index,1)>1

    disp('Calculating dihedrals')

    % Ensure bond_list is sorted with smaller index first
    Bond_index(:,1:2) = sort(Bond_index(:,1:2), 2);

    % Get the total number of atoms
    N_atom = max(Bond_index(:));

    % Initialize adjacency list
    adjacency_list = cell(N_atom, 1);

    % Build adjacency list from bond information
    for i = 1:size(Bond_index, 1)
        atom1 = Bond_index(i, 1);
        atom2 = Bond_index(i, 2);

        adjacency_list{atom1} = [adjacency_list{atom1}, atom2];
        adjacency_list{atom2} = [adjacency_list{atom2}, atom1];
    end

    % Initialize dihedral list
    Dihedral_index = [];

    % Loop over all bonds to identify dihedrals
    for i = 1:size(Bond_index, 1)
        atom2 = Bond_index(i, 1);
        atom3 = Bond_index(i, 2);

        % Neighbors of atom2 excluding atom3
        neighbors2 = adjacency_list{atom2};
        neighbors2(neighbors2 == atom3) = [];

        % Neighbors of atom3 excluding atom2
        neighbors3 = adjacency_list{atom3};
        neighbors3(neighbors3 == atom2) = [];

        % Loop over possible atom1 and atom4 to form dihedrals
        for atom1 = neighbors2
            for atom4 = neighbors3
                % Form the dihedral
                dihedral = [atom1, atom2, atom3, atom4];

                % Enforce a consistent ordering to avoid duplicates
                % (e.g., atom1 < atom4)
                if atom1 < atom4
                    Dihedral_index = [Dihedral_index; dihedral];
                else
                    Dihedral_index = [Dihedral_index; dihedral([4, 3, 2, 1])];
                end
            end
        end
    end

    nDihedrals=size(Dihedral_index,2);

    if nDihedrals>0
        % Remove duplicate dihedrals
        Dihedral_index = unique(Dihedral_index, 'rows');

        [Y,I]=sort(Dihedral_index(:,3));
        Dihedral_index=Dihedral_index(I,:);
        [Y,I]=sort(Dihedral_index(:,2));
        Dihedral_index=Dihedral_index(I,:);

        Pairlist=Dihedral_index(:,[1 4]);
        Pairlist=unique(Pairlist, 'rows');
        Pairlist(ismember(Pairlist, Bond_index(:,1:2), 'rows'),:)=[];
    else
        Dihedral_index =[];
        Pairlist =[];
    end

end

assignin('caller','dist_matrix',dist_matrix);
assignin('caller','overlap_index',overlap_index);
assignin('caller','Bond_index',Bond_index);
assignin('caller','Angle_index',Angle_index);
assignin('caller','nBonds',nBonds);
assignin('caller','nAngles',nAngles);
assignin('caller','nDihedrals',nDihedrals);
assignin('caller','Dihedral_index',Dihedral_index);
assignin('caller','Pairlist',Pairlist);

end