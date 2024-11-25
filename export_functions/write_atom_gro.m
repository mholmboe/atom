%% write_atom_gro.m
% * This function writes a Gromacs .gro file.
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_gro(atom,Box_dim,filename_out) % Basic input arguments
%
function write_atom_gro(atom,Box_dim,filename_out)

if regexp(filename_out,'.gro') ~= false
    filename_out = filename_out;
else
    filename_out = strcat(filename_out,'.gro');
end

if numel(Box_dim)==1
    Box_dim(1)=Box_dim(1);
    Box_dim(2)=Box_dim(1);
    Box_dim(3)=Box_dim(1);
elseif numel(Box_dim)==6
    Box_dim=Cell2Box_dim(Box_dim);
end

%% Box vectors for the .gro format is (free format, space separated reals), values:
% v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
% the last 6 values may be omitted (they will be set to zero) when all angles are 90
% GROMACS only supports boxes with v1(y)=v1(z)=v2(z)=0.

%% Box matrix
% v1(x) v2(x) v3(x)    v1(x) v2(x) v3(x)
% v1(y) v2(y) v3(y) == 0     v2(y) v3(y)
% v1(z) v2(z) v3(z)    0     0     v3(z)

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
