%% reorder_atom_gro.m
% * This function reorders the atoms in a .gro file
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se


%% Examples
% * reorder_atom_gro(atom,[35	34	33	30	27	24	23	26],Box_dim,'out.gro')

function reorder_atom_gro(atom,atomlist,Box_dim,filename_out)
%

% For TCH and from ATB-1 to ATB-2.1
% atomlist=[35	34	33	30	27	24	23	26	39	29	38	18	16	21	20	19	15	37	28	13	36	12	11	10	8	32	31	22	25	17	14	9	6	7	4	3	2	1	5];

nAtoms=size([atom.x],2);
Atom_section=cell(nAtoms,10);

nResidues=max([atom.molid]);
nResidueatoms=nAtoms/nResidues;

atomlist_full=[];
for i=1:nResidues;
    atomlist_full=[atomlist_full (atomlist+(i-1)*nResidueatoms)];
end

fid = fopen(strcat(filename_out,'.gro'), 'wt');
fprintf(fid, '%s\n','Created in Matlab');
fprintf(fid, '%-5i\n',nAtoms);

if isnan(atom(1).vx) || atom(1).vx == 0
    for i = 1:nAtoms
        Atom_section(1:7) = [atom(i).molid, atom(atomlist_full(i)).resname, atom(atomlist_full(i)).type, atom(i).index, atom(atomlist_full(i)).x/10, atom(atomlist_full(i)).y/10, atom(atomlist_full(i)).z/10];
        fprintf(fid, '%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n', Atom_section{1:7});
    end
else
    % Untested
    for i = 1:nAtoms
        Atom_section(1:10) = [atom(i).molid, atom(atomlist_full(i)).resname, atom(atomlist_full(i)).type, atom(i).index, atom(atomlist_full(i)).x/10, atom(atomlist_full(i)).y/10, atom(atomlist_full(i)).z/10, atom(atomlist_full(i)).vx/10, atom(atomlist_full(i)).vy/10, atom(atomlist_full(i)).vz/10];
        fprintf(fid, '%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n', Atom_section{1:10});
    end
end
if size(Box_dim,2) == 3 || Box_dim(4) == 0 && Box_dim(5) == 0 && Box_dim(6) == 0
    fprintf(fid, '%10.5f%10.5f%10.5f\n',Box_dim(1:3)/10);
else
    fprintf(fid, '%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',Box_dim/10);
end
fprintf(fid, '\n');
fclose(fid);
disp('.gro structure file written')
