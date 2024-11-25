%% composition_atom.m
% * This function scans the composition of the atom struct
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
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
composition.natoms = resname_natoms;

Atom_types=[];Atom_numbers=[];Atom_charge=[];
Atom_label=unique([atom.type]);
for i=1:length(unique([atom.type]))
    Atom_types=[Atom_types Atom_label(i)];
    Atom_numbers=[Atom_numbers sum(ismember([atom.type],Atom_label(i)))];
    if isfield(atom,'charge')
        Atom_charge=[Atom_charge mean([atom(ismember([atom.type],Atom_label(i))).charge])];
    end
end
composition.Atom_types=Atom_types
composition.Atom_numbers=Atom_numbers
if isfield(atom,'charge')
    composition.Atom_ave_charge=Atom_charge
    disp('Sum of all charges')
    round([sum([composition.Atom_ave_charge].*[composition.Atom_numbers])     sum([atom.charge])],7)
end

Atom_types=[Atom_types;num2cell(Atom_numbers)];

assignin('caller','Atom_types',Atom_types);
assignin('caller','composition',composition);
