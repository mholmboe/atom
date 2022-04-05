%% check_atom.m
% * This function analyzes various things of the atom struct
% * atom is the atom struct
% * Box_dim is the box dimension vector
% * This function might be outdated...
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = check_atom(atom,Box_dim,1.25,2.25)
%
function atom = check_atom(atom,Box_dim,rmaxshort,rmaxlong)
nAtoms=size(atom,2);

atom=wrap_atom(atom,Box_dim);

%%%%%%%%%% To analyze composition
composition_atom(atom);
Box_dim(1:end)
composition
bond_angle_atom(atom,Box_dim,rmaxshort,rmaxlong,'more');
Atom_labels=unique([atom.type]);

assignin('caller','composition',composition);
assignin('caller','Dist_matrix',Dist_matrix)
assignin('caller','Bond_index',Bond_index);
assignin('caller','Angle_index',Angle_index);
assignin('caller','nBonds',nBonds);
assignin('caller','nAngles',nAngles);

%delete('./#*'); delete('./temp*');
