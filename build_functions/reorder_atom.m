%% reorder_atom.m
% * This function reorders the atoms in an atom struct. Useful for instance
% * for creating united-atom structures from all-atom structures, by simply
% * ignoring the non-polar H's from the original list of atoms, or
% * reordering the atom struct with respect to residue name or atom type
% *
% * In case of reordering the the atom order, neworder is a [1xn] array of
% * n index numbers with a new order. Else if varargin is 'resname' or
% * 'atomtype', neworder is a cell list of 'stringnames'
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = reorder_atom(atom,[1 3 4 5 6 7 8 10]) % Orders according to index
% # atom = reorder_atom(atom,{'MMT' 'SOL' 'Na'},'resname') % Orders according to resname
% # atom = reorder_atom(atom,{'Na' 'Ow' 'Hw'},'atomtype') % Orders according to atomtype
%
function atom = reorder_atom(atom,neworder,varargin)

if nargin>2
    type=varargin{1};
end

if iscell(neworder)==0
    disp('Will re-order the order of the atoms');
    nResidues=max([atom.molid]) % How many residues do we have?
    nResidueatoms_prev=sum([atom.molid]==1); % How many atoms per residue did we have initially?
    nResidueatoms_current=length(neworder); % How many atoms per residue did will we output?
    nAtoms=length(neworder)*nResidues; % How many atoms in total will we output? Can be lower than the initial number

    ordered_atom=atom(1:nAtoms); % Create an atom struct to overwrite
    for i=1:nResidues
        ordered_atom(1+(i-1)*nResidueatoms_current:i*nResidueatoms_current)=atom(neworder+(i-1)*nResidueatoms_prev);
    end

    atom=update_atom(ordered_atom);

    %     % Write the new structure to a .gro file
    %     write_atom_gro(newatom,Box_dim,'reordered.gro')

elseif iscell(neworder)==1 && strncmpi(type,'resname',3)

    resnames=[];ordered_atom=[];
    if numel(neworder)~=numel(unique([atom.resname]))
        disp('You have not issued a complete list of resnames!');
        disp('only this resname was found...');
        unique([atom.resname])
        pause(1)
    else
        if sum(ones(1,numel(neworder))-strcmp(sort(unique(upper(neworder))),sort(unique(upper([atom.resname])))))~=0
            disp('Make sure you have the resnames spelled right..');
        end
    end

    for i=1:numel(neworder)
        resnames(i).ind = find(ismember([atom.resname],neworder(i)));
    end
    for i=1:numel(neworder)
        ordered_atom=[ordered_atom atom(resnames(i).ind)];
    end
    atom=update_atom(ordered_atom);

elseif iscell(neworder)==1 && strncmpi(type,'atomtype',6) || strncmpi(type,'type',4)

    atomtype=[];ordered_atom=[];
    if numel(neworder)~=numel(unique([atom.type]))
        disp('You have not issued a list of atomtypes that fully match the atom struct!');
        if numel(setdiff(unique([atom.type]),unique(neworder)))>0
            setdiff(unique([atom.type]),unique(neworder))
            unique([atom.type])
            pause
        end
    else
        if sum(ones(1,numel(neworder))-strcmp(sort(unique(upper(neworder))),sort(unique(upper([atom.type])))))~=0
            disp('Make sure you have the atomtypes spelled correclty..');
        end
    end

    % neworder(ismember(neworder,{'H' 'O'}))=[];
    % neworder(end+1)={'O'};neworder(end+1)={'H'};

    for i=1:numel(neworder)
        atomtype(i).ind = find(ismember([atom.type],neworder(i)));
    end
    for i=1:numel(neworder)
        ordered_atom=[ordered_atom atom(atomtype(i).ind)];
    end
    atom=update_atom(ordered_atom);

end

assignin('caller','neworder',neworder)

end


