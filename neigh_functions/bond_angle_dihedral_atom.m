%% bond_angle_dihedral_atom.m
% * This function tries to find all bonds, angles and the dihedral angles
% * of the atom struct.
% * Box_dim is the box dimension vector
%
%% Version
% 2.11
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

Dihedral_index=[];
if size(Angle_index,1)>1
    
    disp('Calculating dihedrals')
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
                    Dihedral_index(d,1:5)=[Ax2(i,1) Ax2(i,2) Ax2(i,3) Ax2(j,3) round(theta,2)];
                else
                    Dihedral_index(d,1:5)=[Ax2(j,3) Ax2(i,3) Ax2(i,2) Ax2(i,1) round(theta,2)];
                end
                d=d+1;
            end
        end
        if mod(i,1000)==1
            if i-1>0
                i-1
            end
        end
    end
    
end

nDihedrals=size(Dihedral_index,2);

if nDihedrals>0
    [Y,I] = sort(Dihedral_index(:,2));
    Dihedral_index = Dihedral_index(I,:);
    Dihedral_index = unique(Dihedral_index,'rows','stable');
    Dihedral_index(~any(Dihedral_index,2),:) = [];
else
    Dihedral_index =[];
end

try
    assignin('caller','Ax2',Ax2);
    assignin('caller','dist_matrix',dist_matrix);
    assignin('caller','overlap_index',overlap_index);
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','Angle_index',Angle_index);
    assignin('caller','Dihedral_index',Dihedral_index);
    assignin('caller','nBonds',nBonds);
    assignin('caller','nAngles',nAngles);
    assignin('caller','nDihedrals',nDihedrals);
catch
    assignin('caller','dist_matrix',dist_matrix);
    assignin('caller','overlap_index',overlap_index);
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','Angle_index',Angle_index);
    assignin('caller','Dihedral_index',Dihedral_index);
    assignin('caller','nBonds',nBonds);
    assignin('caller','nAngles',nAngles);
    assignin('caller','nDihedrals',nDihedrals);
end

