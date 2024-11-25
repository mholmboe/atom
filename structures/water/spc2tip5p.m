%% spc2tip5p.m
% * This function converts a .gro or .pdb file with spc water to some tip5p
% water
% * The coordinate of the new MW center is set to OW coordinates, thus must
% be properly energy minimized this could easliy be avoided if using an 
% unwrapped structure...
% * The coordinate of the new MW centers is set to OW coordinates, thus 
% must be properly energy minimized this could easliy be avoided if using 
% an unwrapped structure...
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Dependent functions
% * import_atom.m
% * write_atom_gro.m
% * update_atom.m
% 
%% Examples
% * tip5p_atom = spc2tip5p('SOL.pdb')
% * tip5p_atom = spc2tip5p('SOL.gro')

function tip5p_atom = spc2tip5p(filename)

atom=import_atom_gro(filename); % Import the gro file
nAtoms=size(atom,2); % Total number of atoms
Ow_ind=find(strcmpi([atom.type],'Ow')); % find the index of the water-Oxygens
SOL_ind=sort([Ow_ind Ow_ind+1 Ow_ind+2]); % find the index of the water-Hydrogens
SOL_atom=atom(SOL_ind); % Extract the water into a new atom struct (why not just use atom.resname?)
[SOL_atom.resname]=deal({'SOL'}); % Set the new resname
nSOL=size(SOL_atom,2); % Number of water atoms

new_SOL=[SOL_atom(1:3:end) SOL_atom(1:3:end) SOL_atom(1:3:end) SOL_atom(1:3:end) SOL_atom(1:3:end)]; % Create a new struct
new_SOL(1:5:end)=SOL_atom(1:3:end); % assign Ow -> Ow
new_SOL(2:5:end)=SOL_atom(2:3:end); % assign Hw -> Hw
new_SOL(3:5:end)=SOL_atom(3:3:end); % assign Hw -> Hw
new_SOL(4:5:end)=SOL_atom(1:3:end); % assign Ow -> Mw
new_SOL(5:5:end)=SOL_atom(1:3:end); % assign Ow -> Mw

[new_SOL(2:5:end).type]=deal({'HW1'}); % rename Ow -> Mw
[new_SOL(3:5:end).type]=deal({'HW2'}); % rename Ow -> Mw
[new_SOL(4:5:end).type]=deal({'LP1'}); % rename Ow -> Mw
[new_SOL(5:5:end).fftype]=deal({'LP2'}); % rename Ow -> Mw

% Put the new_atom struct back again into the original atom struct
if SOL_ind(end) == nAtoms; 
    tip5p_atom=[atom(1:SOL_ind(1)-1) new_SOL];
else
    tip5p_atom=[atom(1:SOL_ind(1)-1) new_SOL atom(SOL_ind(end)+1:nAtoms)];
end

% Update the atoms index and stuff
tip5p_atom=update_atom(tip5p_atom);

% Write the new .gro file
write_atom_gro(tip5p_atom,Box_dim,strcat('tip5p_',filename))

