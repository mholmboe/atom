%% number_type.m
% * This function numbers the atom types, like H1, H2, H3...
% * from the atom struct
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = number_type(atom) % Basic input arguments
% * atom = number_type(atom,1) % 1 for add numbers, O for remove the numbers in the atom type names

function atom = number_type(atom,varargin)

Molids=unique([atom.molid]);

if nargin==1
    atom=element_atom(atom);
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
elseif nargin>1
    if varargin{1}==1
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
    else
        for i=1:numel(Molids)
            Atom_labels=unique([atom(ismember([atom.molid],Molids(i))).type]);
            for j=1:numel(Atom_labels)
                n=1;
                if ~strncmpi(Atom_labels(j),'Ow',2)
                    ind=find([atom.molid]==Molids(i)&ismember([atom.type],Atom_labels(j)));
                    for k=1:numel(ind)
                        atom(ind(k)).type=regexprep([atom(ind(k)).type],'\d+$','');
                        n=n+1;
                    end
                end
            end
        end

    end

end

end



