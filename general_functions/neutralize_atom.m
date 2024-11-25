%% neutralize_atom.m
% * This function appends a 0 to all atomtype names and will also set the 
% charge (if the field charge exist) to zero (0).
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = neutralize_atom(atom)
%
function atom = neutralize_atom(atom)

atom=element_atom(atom);
nAtoms=size(atom,2);

for i=1:nAtoms
    regexprep([atom(i).type],'[\d"]','');
    [atom(i).type]=strcat(atom(i).type{1}(end),{'0'});
end

if isfield(atom,'charge')
    [atom.charge]=deal(0);
end
