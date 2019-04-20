%% write_atom_gro.m
% * This function writes a gro file. Does it also write velocities?
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_gro(atom,Box_dim,filename_out)
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
end

nAtoms=length(atom);
Atom_section=cell(nAtoms,10);
fid = fopen(filename_out, 'wt');
fprintf(fid, '%s\n','Created in Matlab');
fprintf(fid, '%-5i\n',nAtoms);

if sum(find(isnan([atom.vx]))) || atom(1).vx == 0
    for i = 1:nAtoms
        Atom_section(1:7) = [atom(i).molid, atom(i).resname, atom(i).type, atom(i).index, atom(i).x/10, atom(i).y/10, atom(i).z/10];
        fprintf(fid, '%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n', Atom_section{1:7});
    end
else
    % Untested
    for i = 1:nAtoms
        Atom_section(1:10) = [atom(i).molid, atom(i).resname, atom(i).type, atom(i).index, atom(i).x/10, atom(i).y/10, atom(i).z/10, atom(i).vx/10, atom(i).vy/10, atom(i).vz/10];
        fprintf(fid, '%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n', Atom_section{1:10});
    end
end
Box_dim
if size(Box_dim,2) == 3 || Box_dim(6) == 0 && Box_dim(8) == 0 && Box_dim(9) == 0
    fprintf(fid, '%10.5f%10.5f%10.5f\n',Box_dim(1:3)/10);
else
    fprintf(fid, '%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',Box_dim/10);
end
%fprintf(fid, '\n');
fclose(fid);
disp('.gro structure file written')
