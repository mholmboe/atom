%% number_type.m
% * This function numbers the atom types, like H1, H2, H3...
% * from the atom struct
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% * number_type(atom)


function atom = number_type(atom,varargin)

atom=element_atom(atom);

Molids=unique([atom.molid]);

for i=1:numel(Molids)
    Atom_labels=unique([atom(ismember([atom.molid],Molids(i))).type]);
    for j=1:numel(Atom_labels)
        n=1;
        if ~strncmpi(Atom_labels(j),'Ow',2)
            ind=find([atom.molid]==Molids(i)&ismember([atom.type],Atom_labels(j)));
            for k=1:numel(ind)
                atom(ind(k)).type=strcat([atom(ind(k)).type],num2str(n));
                n=n+1;
            end
        end
    end
end





