%% noupdate_atom.m
% * This function updates atom index but not the the molid index in the 
% atom struct
% * Multiple atom structs can be also concatenated by using this format
% atom = update_atom({atom1 atom2 atom3})
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = noupdate_atom(atom) % Basic input arguments
% # atom = noupdate_atom({atom1 atom2 atom3}) % Will append atom1 and atom2 and atom3 and update their MolID's, respectively
%
function atom = noupdate_atom(atom,varargin)

if iscell(atom)
    
    % In case first struct is empty
    if size(atom,2)>1 && size(atom{1},2)==0 % Was 1?
        for i=2:size(atom,2)
            newatom{i-1}=atom{i};
        end
        atom=newatom;
    end
    
    size_ind=zeros(size(atom,2),1);
    for i=1:size(atom,2)
        if size(atom{i},2)>0
            size_ind(i)=1;
        end
    end
    
    keepfieldnames=fieldnames(atom{1}); % Orig line
    if size(atom,2) > 1
        for i=1:size(atom,2)
            if numel(atom{i})
                keepfieldnames=intersect(keepfieldnames,fieldnames(atom{i}));
            end
        end
        for i=1:size(atom,2)
            if numel(atom{i})
                rmfieldnames=setdiff(fieldnames(atom{i}),keepfieldnames);
                for j=1:numel(rmfieldnames)
                    atom{i}=rmfield(atom{i},rmfieldnames(j));
                end
            end
        end
    end
    
    atom=atom(logical(size_ind));
    atom_Tot=atom{1}; % Orig line
    if size(atom,2) > 1
        for i=2:size(atom,2)
            atom_temp=atom{i};
            if size(atom_temp,2)>0
                if numel(unique([atom_temp.molid]))==1
                    molid=num2cell(([atom_temp.molid]+[atom_Tot(end).molid]));
                    [atom_temp.molid]=deal(molid{:});
                elseif numel(unique([atom_temp.molid]))==size(atom_temp,2)
                    molid=num2cell(([1:size(atom_temp,2)]+[atom_Tot(end).molid]));
                    [atom_temp.molid]=deal(molid{:});
                end
            end
            atom_Tot=[atom_Tot atom_temp];
        end
    end
    atom=atom_Tot;
end

nAtoms=size(atom,2);

atom=order_attributes(atom); % Order all attributes

% assignin('base','atom1',atom);

% atom = resname_atom(atom);

assignin('caller','nAtoms',nAtoms);

end

