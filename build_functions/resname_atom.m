%% resname_atom.m
% * This function tries to guess the resname based in the atom types
%
%% Version
% 2.081
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
ion={'Li','LI','Li+','LI+','Na','NA','Na+','NA+','K','K+','Rb','RB','Rb+','RB+','Cs','CS','Cs+','CS+',...
     'Ca','CA','Ca2+','CA2+','Cu','CU','Cu2+','CU2+','Ni','NI','Ni2+','NI2+',...
     'Zn','ZN','Zn2+','ZN2+','Sr','SR','Sr2+','SR2+','Ba','BA','Ba2+','BA2+','F','F-','Cl','CL','Cl-','CL-',...
     'Br','BR','Br-','BR-','I','I-'}; % 'Mg','MG','Mg2+','MG2+',
ION=upper(ion);ion=[ion ION];

Sol_ind=sort([find(strncmpi([atom.type],sol(1),2)) find(strncmpi([atom.type],sol(2),2))]);
[atom(Sol_ind).resname]=deal({'SOL'});
Ion_ind=find(ismember([atom.type],ion));
[atom(Ion_ind).resname]=deal({'ION'});
% noSol_ind=setdiff([atom.index],Sol_ind);
