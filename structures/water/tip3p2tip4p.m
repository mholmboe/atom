%% tip3p2tip4p.m
% * This function converts a .gro or .pdb file with spc water to some tip4p
% water
% * The coordinate of the new MW center is set to OW coordinates, thus must
% be properly energy minimized this could easliy be avoided if using an 
% unwrapped structure...
% * The coordinate of the new MW centers is set to OW coordinates, thus 
% must be properly energy minimized this could easliy be avoided if using 
% an unwrapped structure...
%
%% Version
% 2.082
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
% * tip4p_atom = tip3p2tip4p('SOL.pdb')
% * tip4p_atom = tip3p2tip4p('SOL.gro')
% * tip4p_atom = tip3p2tip4p('SOL.pdb')
% * tip4p_atom = tip3p2tip4p('SOL.gro')

function tip4p_atom = tip3p2tip4p(filename)

atom=import_atom(filename); % Import the .pdb/.gro file
nAtoms=size(atom,2); % Total number of atoms
Ow_ind=find(strcmpi([atom.type],'Ow')); % find the index of the water-Oxygens
SOL_ind=sort([Ow_ind Ow_ind+1 Ow_ind+2]); % find the index of the water-Hydrogens
SOL_atom=atom(SOL_ind); % Extract the water into a new atom struct (why not just use atom.resname?)
[SOL_atom.resname]=deal({'SOL'}); % Set the new resname
nSOL=size(SOL_atom,2); % Number of water atoms

new_SOL=[SOL_atom(1:3:end) SOL_atom(1:3:end) SOL_atom(1:3:end) SOL_atom(1:3:end)]; % Create a new struct
new_SOL(1:4:end)=SOL_atom(1:3:end); % assign Ow -> Ow
new_SOL(2:4:end)=SOL_atom(2:3:end); % assign Hw -> Hw
new_SOL(3:4:end)=SOL_atom(3:3:end); % assign Hw -> Hw
new_SOL(4:4:end)=SOL_atom(1:3:end); % assign Ow -> Mw

[new_SOL(4:4:end).type]=deal({'MW'}); % rename Ow -> Mw
[new_SOL(4:4:end).fftype]=deal({'MW'}); % rename Ow -> Mw

% Put the new_atom struct back again into the original atom struct
if SOL_ind(end) == nAtoms; 
    tip4p_atom=[atom(1:SOL_ind(1)-1) new_SOL];
else
    tip4p_atom=[atom(1:SOL_ind(1)-1) new_SOL atom(SOL_ind(end)+1:nAtoms)];
end

% Update the atoms index and stuff
tip4p_atom=update_atom(tip4p_atom);

% Write the new .gro file
write_atom_gro(tip4p_atom,Box_dim,strcat('tip4p_',filename))

