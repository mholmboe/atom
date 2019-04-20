%% update_atom.m
% * This function updates the molid index and the atoms index in the atom struct
% * Multiple atom structs can be also concatenated by using this format atom = update_atom({atom1 atom2 atom3})
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% #
% # atom = update_atom(atom)
% # atom = update_atom({atom1 atom2 atom3})
%
function atom = update_atom(atom)

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
            keepfieldnames=intersect(keepfieldnames,fieldnames(atom{i}));
        end
        for i=1:size(atom,2)
            rmfieldnames=setdiff(fieldnames(atom{i}),keepfieldnames);
            for j=1:numel(rmfieldnames)
                atom{i}=rmfield(atom{i},rmfieldnames(j))
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

nAtoms=size([atom.x],2);

MolID=[atom.molid];

if length(MolID) < nAtoms
    MolID=[MolID;[MolID(end)+1:nAtoms]'];
end

% Update Molid and index
nmol=1;first_in=[1];last_in=[];
for i=1:nAtoms
    if i > 1 && MolID(i) ~= MolID(i-1) | strcmp([atom(i).resname],[atom(i-1).resname]) == 0
        nmol=nmol+1;
        atom(i).molid=nmol;
        first_in(atom(i).molid,1)=i; last_in(atom(i).molid-1,1)=i-1;
    elseif i > 1
        atom(i).molid=atom(i-1).molid;
    elseif i == 1
        atom(i).molid=1;
    end
    atom(i).index=mod(i,100000);
end

if numel(fieldnames(atom))~=10
    defaultAttributes={'molid' 'resname' 'type' 'fftype' 'index' 'neigh' 'bond' 'angle' 'x' 'y' 'z' 'vx' 'vy' 'vz' 'xfrac' 'yfrac' 'zfrac' 'element' 'mass' 'Mw' 'COM_x' 'COM_y' 'COM_z' 'charge' 'bv' 'mean_bv' 'valence' 'Rdiff' 'atnum'};
    atomAttributes=fieldnames(atom)';
    indDefault=find(ismember(defaultAttributes,atomAttributes));
    defaultAttributes=defaultAttributes(indDefault);
    ind_atom=find(ismember(atomAttributes,defaultAttributes));
    atomAttributes=atomAttributes(ind_atom);
    atom=orderfields(atom,unique({defaultAttributes{:} atomAttributes{:}},'stable'));
end
% assignin('base','atom1',atom);

% atom = resname_atom(atom);

assignin('caller','nAtoms',nAtoms);

