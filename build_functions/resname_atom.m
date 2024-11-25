%% resname_atom.m
% * This function tries to guess the resname based in the atom types
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = resname_atom(atom) % Basic inout arguments
%
function atom = resname_atom(atom,varargin)

nAtoms=size(atom,2);
Atom_label=sort(unique([atom.type]));

% clay = {'H','Oh','O','Op','Ob','Omg', 'Oalt','Odsub','Ohmg','Oalh','Oalhh','Oalsi','Osih','Si','Al','Alt','Mgo','Mgh','Cao','Cah','Feo','Lio'};
% From Cygan, 2004  = {'h*','ho','o*','oh','ob','obos','obts', 'obss', 'ohs','OAlH','OAlH2','OAlSi','OSiH','st','ao','Alt','mgo','mgh','cao','cah','feo','lio','Li','Na','K','Cs','Mg','Ca','Sr','Ba','Cl','Br'};

sol={'Hw','Ow'};
ion={'Li','Li+','LI+','Na','NA','Na+','NA+','K','K+','Rb','RB','Rb+','RB+','Cs','CS','Cs+','CS+',...
    'Ca','CA','Ca2+','CA2+','Cu','CU','Cu2+','CU2+','Ni','NI','Ni2+','NI2+',...
    'Zn','ZN','Zn2+','ZN2+','Sr','SR','Sr2+','SR2+','Ba','BA','Ba2+','BA2+','F-','Cl','CL','Cl-','CL-',...
    'Br','BR','Br-','BR-','I','I-'}; % 'Mg','MG','Mg2+','MG2+',
ION=upper(ion);ion=[ion ION];

SOL_ind=sort([find(strncmpi([atom.type],sol(1),2)) find(strncmpi([atom.type],sol(2),2))]);
[atom(SOL_ind).resname]=deal({'SOL'});
ION_ind=find(ismember([atom.type],ion));
[atom(ION_ind).resname]=deal({'ION'});
% noSol_ind=setdiff([atom.index],Sol_ind);

index=num2cell(1:size(atom,2));
[atom.index]=deal(index{:});
ind=~ismember([atom.index],[SOL_ind ION_ind]);
if numel(ind)>0
    if numel([atom(ind(1)).resname])<2
        [atom(~ismember([atom.index],[SOL_ind ION_ind])).resname]=deal({'MIN'});
    end
end

if nargin>1
    resname=varargin{1};
    if ~iscell(resname)
        resname={resname};
    end
    [atom(~ismember([atom.index],[SOL_ind ION_ind])).resname]=deal(resname);
end

end