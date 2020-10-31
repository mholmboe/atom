%% composition_atom.m
% * This function scans the composition of the atom struct
%
%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = composition_atom(atom)

function atom = composition_atom(atom)

resname_labels=unique([atom.resname]);
resname_natoms=zeros(1,length(resname_labels));resname_nresidues=zeros(1,length(resname_labels));
for i=1:length(resname_labels)
    resname_ind=ismember([atom.resname],resname_labels(i));
    resname_natoms(i)=sum(resname_ind);
    resname_nresidues(i)=length(unique([atom(resname_ind).molid]));
end

composition.resnames = resname_labels;
composition.nresidues = resname_nresidues;
composition.natoms = resname_natoms

Atom_types=[];Atom_numbers=[];
Atom_label=unique([atom.type]);
for i=1:length(unique([atom.type]))
    Atom_types=[Atom_types Atom_label(i)];
    Atom_numbers=[Atom_numbers sum(ismember([atom.type],Atom_label(i)))];
end
Atom_types
Atom_numbers

Atom_types=[Atom_types;num2cell(Atom_numbers)];

assignin('caller','Atom_types',Atom_types);
assignin('caller','composition',composition);
