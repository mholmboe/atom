%% resname_atom.m
% * This function tries to guess the resname based in the atom types
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = resname_atom(atom)
%
function atom = resname_atom(atom)

nAtoms=size(atom,2);
Atom_label=sort(unique([atom.type]));

% clay = {'H','Oh','O','Op','Ob','Omg', 'Oalt','Odsub','Ohmg','Oalh','Oalhh','Oalsi','Osih','Si','Al','Alt','Mgo','Mgh','Cao','Cah','Feo','Lio'};
% From Cygan, 2004  = {'h*','ho','o*','oh','ob','obos','obts', 'obss', 'ohs','OAlH','OAlH2','OAlSi','OSiH','st','ao','Alt','mgo','mgh','cao','cah','feo','lio','Li','Na','K','Cs','Mg','Ca','Sr','Ba','Cl','Br'};

sol={'Hw','Ow'};
ion={'Li','Na','Na+','K','K+','Rb','Cs','Mg','Ca','Cu','Ni','Zn','Sr','Ba','F','Cl','Br','I'};
ION=upper(ion);ion=[ion ION];

Sol_ind=sort([find(strncmpi([atom.type],sol(1),2)) find(strncmpi([atom.type],sol(2),2))]);
[atom(Sol_ind).resname]=deal({'SOL'});
Ion_ind=find(ismember([atom.type],ion));
[atom(Ion_ind).resname]=deal({'ION'});
% noSol_ind=setdiff([atom.index],Sol_ind);
