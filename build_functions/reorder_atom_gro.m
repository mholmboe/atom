%% reorder_atom_gro.m
% * This function reorders the atoms in a .gro file. Useful for creating
% * united-atom structures from all-atom structures. The script assumes there
% * is only one type of molecule in the atom struct.
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # reorder_atom_gro(atom,[35	34	33	30	27	24	23	26],Box_dim,'out.gro')
% # reorder_atom_gro(atom,[27 26 23 20 17 14 11 8 5 1],Box_dim,'ua_oct.gro')

function reorder_atom_gro(atom,atomlist,Box_dim,filename_out)

nResidues=max([atom.molid]); % How many residues are we talking about?
nAtoms=length(atomlist)*nResidues; % How many atoms in total will we output?
Atom_section=cell(nAtoms,10); 
nResidueatoms_prev=sum([atom.molid]==1); % How many atoms per residue did we have initially?
nResidueatoms_current=length(atomlist); % How many atoms per residue did will we output?

atomlist_full=[];
newatom=atom(1:nAtoms); % Create an atom struct to overwrite
for i=1:nResidues
    atomlist_full=[atomlist_full (atomlist+(i-1)*nResidueatoms_prev)];
    newatom(1+(i-1)*nResidueatoms_current:i*nResidueatoms_current)=atom(atomlist+(i-1)*nResidueatoms_prev);
end

% Write the new structure to a .gro file
write_atom_gro(newatom,Box_dim,filename_out)
