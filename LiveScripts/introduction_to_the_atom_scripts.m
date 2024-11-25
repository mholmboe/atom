%% Creating and manipulating atoms/molecules/structures in Matlab using the atom scripts
% Molecular simulations need physically sound starting structure/configurations 
% (especially of multi-component systems), as well as a correct molecular and 
% forcefield specific molecular topology (information of the atomic, bonding and 
% angle order, atom type info and forcefield parameters etc..)
% 
% 
% 
% Gromacs own tools can help in many situations, especially with protein or 
% lipid based systems, see…
% 
% |gmx  {editconf, solvate, insert-molecules, genion, pdb2gmx, x2top and so 
% on}|
% 
% 
% 
% VMD with various plugins – the psfgen, solvate, topotools or inorganic builder 
% plugins are some useful alternatives. GUI and tk-console versions often avail.
% 
% 
% 
% However… since no code will always work for every type of system/simulation 
% setup, being able to write one’s own ’tools’ is sometimes neccessary
% 
% 
%% When and why?
% To perform isomorphic subsitutions
% 
% To find bonds with variable cutoff’s
% 
% Automatically find bonds across the PBC (can VMD’s topotool do this?)
% 
% To find and assign the correct atomtypes for Clayff, Interface
% 
% etc..
% 
% For the purpose of automizing/enabling efficient construction of complex and 
% anisotropic (geo)chemical systems, and the corresponding topological info, we 
% use Matlab and matlab struct variable.
% 
% 
%% Why?
% Perform isomorphic subsitutions
% 
% Find bonds with variable cutoff’s
% 
% Find bonds across the PBC
% 
% Find and assign the correct atomtypes for Clayff, Interface etc
% 
% Easily build multicomponent systems {mineral/ions/water/organics}
% 
% Manipulate the system on the atomic, molecule or molecular type level
% 
% 
%% How?
% Examples of the atom functions! Using the struct variable ’atom’ in matlab

atom = import_atom(filename) % xyz/pdb/gro
atom = write_atom(arguments) % xyz/pdb/gro/mol2/pqr and lammps/gmx topologies
atom = solvate_atom(arguments)
atom = merge_atom(arguments)
atom = translate_atom(arguments)
atom = COM_atom(arguments)
%% 
% +50 other functions…
% 
% 
%% *Abouth paths…*
% The VMD path on your computer. Useful to call VMD from within Matlab

function PATH2VMD = PATH2VMD()
    % Add your own path to VMD (this works on my Mac..)
    PATH2VMD = '/Applications/VMD\ 1.9.2.app/Contents/MacOS/startup.command';
end
%% 
% 
% 
% 
% 
% 
% 
%