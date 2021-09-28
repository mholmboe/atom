%% dist_matrix_noPBC_atom.m
% * This function calculates the distance matrix from the atom struct, or 
% the distances between atom1 and atom2
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # dist_matrix = dist_matrix_noPBC_atom(atom1) % Basic input arguments
% # dist_matrix = dist_matrix_noPBC_atom(atom1,atom2) % Calculates the distance matrix between sites in atom1 and in atom2
%
function dist_matrix = dist_matrix_noPBC_atom(atom1,varargin) % ,atom2,Box_dim); % or % ,Box_dim);

format compact;

if nargin==2
    atom2=varargin{1};
else
    atom2=atom1;
end

nAtoms1=size(atom1,2);
nAtoms2=size(atom2,2);

XYZ1=single([[atom1.x]' [atom1.y]' [atom1.z]']);
XYZ2=single([[atom2.x]' [atom2.y]' [atom2.z]']);

dist_matrix = single(zeros(nAtoms2,nAtoms1)); % use of single instead of double 
X_dist = dist_matrix;
Y_dist = dist_matrix;
Z_dist = dist_matrix;

    for i = 1:size(XYZ1,1)
        %Calculate Distance Components
        rx = XYZ1(i,1) - XYZ2(:,1);

        ry = XYZ1(i,2) - XYZ2(:,2);

        rz = XYZ1(i,3) - XYZ2(:,3);

        r = sqrt( rx(:,1).^2 + ry(:,1).^2 + rz(:,1).^2 ); % distance calc.
        dist_matrix(:,i)=r;
        X_dist(:,i)=rx;
        Y_dist(:,i)=ry;
        Z_dist(:,i)=rz;
        
        if mod(i,1000)==1
            if i > 1
                i-1
            end
        end
    end

% New transposed output
dist_matrix=dist_matrix';
try
    assignin('caller','X_dist',(X_dist)');
    assignin('caller','Y_dist',(Y_dist)');
    assignin('caller','Z_dist',(Z_dist)');
catch
    assignin('base','X_dist',(X_dist)');
    assignin('base','Y_dist',(Y_dist)');
    assignin('base','Z_dist',(Z_dist)');
end

