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
    atom=bond_atom(atom,Box_dim,rmaxlong);
elseif nargin>4
    Box_dim=varargin{1};
    rmaxshort=varargin{2};
    rmaxlong=varargin{3};
    atom=bond_angle_atom(atom,Box_dim,rmaxshort,rmaxlong,'more');
end

Dihedral_index=[];Pairlist=[];Improper_dihedral_index = [];
nDihedrals=size(Dihedral_index,2);
if size(Bond_index,1)>1

    disp('Calculating dihedrals and a pairlist')

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

    if size(Dihedral_index,1)>1
        %    Ax2=[[Angle_index(:,3) Angle_index(:,2) Angle_index(:,1) Angle_index(:,4) Angle_index(:,8:10) Angle_index(:,5:7)]; Angle_index];
        Ax2=[Angle_index(:,[3 2 1 4 8 9 10 5 6 7]); Angle_index];
        d=1;
        for i=1:size(Ax2,1)
            for j=i:size(Ax2,1)
                if isequal([Ax2(i,2) Ax2(i,3)],[Ax2(j,1) Ax2(j,2)])
                    A=cross([Ax2(i,5) Ax2(i,6) Ax2(i,7)],[Ax2(i,8) Ax2(i,9) Ax2(i,10)]);
                    B=cross([Ax2(j,5) Ax2(j,6) Ax2(j,7)],[Ax2(j,8) Ax2(j,9) Ax2(j,10)]);
                    normA=sqrt(sum(A.*A,2));
                    normB=sqrt(sum(B.*B,2));
                    theta=rad2deg(acos(dot(A,B)./(normA*normB)));
                    if Ax2(i,2)<Ax2(i,3)
                        Dihedral_index(d,1:5)=real([Ax2(i,1) Ax2(i,2) Ax2(i,3) Ax2(j,3) round2dec(theta,2)]);
                    else
                        Dihedral_index(d,1:5)=real([Ax2(j,3) Ax2(i,3) Ax2(i,2) Ax2(i,1) round2dec(theta,2)]);
                    end
                    d=d+1;
                end
            end
        end
    end

    disp('Calculating improper dihedrals with planarity check...');

    % Define a threshold for planarity.
    % If abs(planarity) is close to 1 within this tolerance, we treat it as planar.
    planarityThreshold = 1e-1;
    indH=find(strncmp([atom.type],'H',1));
    % adjacency_list_noH=adjacency_list;
    % % for i = 1:size(adjacency_list_noH,1)
    % %     adjacency_list_noH{i}(ismember(adjacency_list_noH{i}, indH)) = [];
    % % end

    % Loop over all atoms to find those with exactly three neighbors
    for central_atom = 1:size(adjacency_list, 1)

        ind3=find(ismember(Angle_index(:,2),central_atom));

        if numel(ind3) == 3 && sum(ismember(atom(central_atom).type,{'CG2O1', 'CG2O3'})) > 0

            % --------------------------------------------------------------
            % 1) Extract coordinates of the central atom and its neighbors
            % --------------------------------------------------------------
            pos_vec=[Angle_index(ind3(1),5:7);...
                Angle_index(ind3(1),8:10);...
                Angle_index(ind3(2),5:7);...
                Angle_index(ind3(2),8:10);...
                Angle_index(ind3(3),5:7);...
                Angle_index(ind3(3),8:10)];
            pos_vec=unique(pos_vec,'rows');

            vec1=pos_vec(1,:);
            vec2=pos_vec(2,:);
            vec3=pos_vec(3,:);

            % --------------------------------------------------------------
            % 2) Compute normals to check planarity
            % --------------------------------------------------------------
            normal1 = cross(vec1, vec2);
            normal2 = cross(vec1, vec3);

            % Guard against very small normals to avoid division by zero
            if norm(normal1) < 1e-12 || norm(normal2) < 1e-12
                continue;  % Degenerate or collinear vectors; skip
            end

            % Planarity measure = cos(angle) between the two normals
            planarity = dot(normal1, normal2) / (norm(normal1) * norm(normal2));

            % If near ±1, the four atoms (center + 3 neighbors) are planar
            if abs(abs(planarity) - 1) < planarityThreshold
                dihedralRow = [central_atom, unique([Angle_index(ind3,[1 3])])'];
                Improper_dihedral_index = [Improper_dihedral_index; dihedralRow];
            end
         end
    end


    % Remove any duplicates and sort the final improper dihedral list
    if ~isempty(Improper_dihedral_index)
        Improper_dihedral_index = unique(Improper_dihedral_index, 'rows');
        Improper_dihedral_index = sortrows(Improper_dihedral_index);
    end

    nImpropers=size(Improper_dihedral_index,2);
end

assignin('caller','dist_matrix',dist_matrix);
%assignin('caller','overlap_index',overlap_index);
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