%% analyze_atom.m
% * This function analyzes various things of the atom struct
% * atom is the atom struct
% * Box_dim is the box dimension vector
% * This function might be outdated...
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom = analyze_atom(atom,Box_dim,1.25,2.25)


function atom = analyze_atom(atom,Box_dim,max_H_dist,max_dist)
nAtoms=size(atom,2);

atom=wrap_atom(atom,Box_dim);

%%%%%%%%%% To analyze composition
composition_atom(atom);
Box_dim(1:3)
composition
bond_angle_atom(atom,Box_dim,max_H_dist,max_dist,'more');
Atom_labels=unique([atom.type]);

assignin('caller','composition',composition);
assignin('caller','Dist_matrix',Dist_matrix)
assignin('caller','Bond_index',Bond_index);
assignin('caller','Angle_index',Angle_index);
assignin('caller','nBonds',nBonds);
assignin('caller','nAngles',nAngles);

%delete('./#*'); delete('./temp*');
