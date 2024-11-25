%% write_atom_dodecahedron_gro.m
% * This function writes a gro file.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_gro(atom,Box_dim,d,type,filename_out) % Basic input
% arguments. d is the corresponding cubic Box size, type is either 'square'
% or 'hex'. See http://manual.gromacs.org/current/reference-manual/algorithms/periodic-boundary-conditions.html
%
function [atom,Box_dim,Cell] = write_atom_dodecahedron_gro(atom,Box_dim,d,type,varargin)

if nargin>4
    filename_out = varargin{1};
    if regexp(filename_out,'.gro') ~= false
        filename_out = filename_out;
    else
        filename_out = strcat(filename_out,'.gro');
    end
end


if numel(Box_dim)==1
    Box_dim(1)=Box_dim(1);
    Box_dim(2)=Box_dim(1);
    Box_dim(3)=Box_dim(1);
end

%% Box vectors for the .gro format is (free format, space separated reals), values:
% v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
% the last 6 values may be omitted (they will be set to zero) when all angles are 90
% GROMACS only supports boxes with v1(y)=v1(z)=v2(z)=0.

%% Box matrix
% v1(x) v2(x) v3(x)    v1(x) v2(x) v3(x)
% v1(y) v2(y) v3(y) == 0     v2(y) v3(y)
% v1(z) v2(z) v3(z)    0     0     v3(z)

if strncmpi(type,'hex',3)
    type='hexagon along xy plane'
    %% Rhombic dodcahdron (xy-hexagon)
    Box_matrix=[...
        d    1/2*d       1/2*d;...
        0    1/2*3^0.5*d 1/6*3^0.5*d;...
        0    0           1/3*6^0.5*d;...
        ];
else
    %% Rhombic dodcahdron (xy-square, gromacs default)
    type='square along xy plane'
    Box_matrix=[...
        d    0   1/2*d;...
        0    d   1/2*d;...
        0    0   1/2*2^0.5*d;...
        ];
end

disp('New rhombic dodecahedron Box')
type
Box_dim=[Box_matrix(1,1) Box_matrix(2,2) Box_matrix(3,3) 0 0 Box_matrix(1,2) 0 Box_matrix(1,3) Box_matrix(2,3)]
Cell=Box_dim2Cell(Box_dim)

if nargin>4

    nAtoms=size(atom,2);
    Atom_section=cell(nAtoms,10);
    fid = fopen(filename_out, 'wt');
    fprintf(fid, '%s\n','Created in Matlab');
    fprintf(fid, '%-5i\n',nAtoms);

    if sum(find(isnan([atom.vx]))) || atom(1).vx == 0
        i=1;
        while i<nAtoms+1
            Atom_section(1:7) = [atom(i).molid, atom(i).resname, atom(i).type, atom(i).index, atom(i).x/10, atom(i).y/10, atom(i).z/10];
            fprintf(fid, '%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n', Atom_section{1:7});
            i=i+1;
        end
    else
        % To include velocities (if they exist) as well... untested
        i=1;
        while i<nAtoms+1
            Atom_section(1:10) = [atom(i).molid, atom(i).resname, atom(i).type, atom(i).index, atom(i).x/10, atom(i).y/10, atom(i).z/10, atom(i).vx/10, atom(i).vy/10, atom(i).vz/10];
            fprintf(fid, '%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n', Atom_section{1:10});
            i=i+1;
        end
    end
    % Box_dim
    if size(Box_dim,2) == 3 || Box_dim(6) == 0 && Box_dim(8) == 0 && Box_dim(9) == 0
        fprintf(fid, '%10.5f%10.5f%10.5f\n',Box_dim(1:3)/10);
    else
        fprintf(fid, '%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',Box_dim/10);
    end
    %fprintf(fid, '\n');
    fclose(fid);
    disp('.gro structure file written')
end
