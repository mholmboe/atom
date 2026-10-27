%% check_clayff_H2Odens.m
% * Check the approx. water density for a clayff system through a semi-empirical relation...
% * atom is the atom struct
% * Box_dim is the box dimension vector [x y z]
%
%% Version
% 3.0.1
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * H2Odens = check_clayff_H2Odens(atom,Box_dim)
%
%% See also
% # <check_H2Odens.html check_H2Odens>
% # <check_clayff_2004_charge.html check_clayff_2004_charge>
% # <check_clayff_charge.html check_clayff_charge>
% # <charge_clayff_2004_atom.html charge_clayff_2004_atom>
% # <charge_clayff_atom.html charge_clayff_atom>
%
function H2Odens = check_clayff_H2Odens(atom,Box_dim)

nSOL=sum(strncmpi([atom.type],'Ow',2));
nSOLatoms=sum(ismember([atom.resname],{'Water' 'water' 'WAT' 'wat' 'SOL' 'sol'}));
nAtomsperSOL=nSOLatoms/nSOL;
V = Box_dim(1)*Box_dim(2)*Box_dim(3); nSolutes=size(atom,2)-nSOL;
VSOL=V-nSolutes/0.092; H2Odens=nSOL/nAtomsperSOL*18.09/6.022E23*1E24/VSOL;
