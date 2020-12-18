%% insert_close_atom.m
% * This function inserts a whole molecule from a structure file or atom_in
% into a region defined by <limits> with a atom (molecule) structure
% * rotate can be a string like 'random', {'random'}, or be used to set
% some angles like [60 90 60]. varargin can be used to assure that one
% atom type is at least some distance above (in z) some other atom type.
%
%% Similar
% fuse_atom
% protonate_atom
% create_atom.m
%
%% Version
% 2.081
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = insert_close_atom(atom,limits,'rotate',r,maxsol)
% * atom = insert_close_atom(atom,limits,[10 20 30],r,maxsol,solute_atom)
% * atom = insert_close_atom(atom,limits,'rotate',r,maxsol,solute_atom,{'C1' 'N1'},0.3)
% * atom = insert_close_atom(atom,limits,'rotate',r,maxsol,solute_atom,[1 4],0.3)

function atom = insert_close_atom(atom_in,limits,rotate,r,nmax,varargin)

atom = insert_close_atom(atom,limits,[10 20 30],r,maxsol,solute_atom)

i=0
while i<MIN_charge
    OH = insert_close_atom(OH,Full_Box_dim,'rotate',3,MIN_charge,System);
    d=dist_matrix_atom(OH,System,Full_Box_dim);
    if min(d1,:)<5
        System = update_atom({System OH});
    end
    i=i+1;
end


