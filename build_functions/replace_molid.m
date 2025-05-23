%% replace_molid.m
% * This function replaces a molecule (by its MolID) in an atom struct with
% a new (single MolID) atom struct by placing the COM of the latter in the
% place of the former
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = replace_molid(in_atom,old_atom,3)

function atom=replace_molid(new_atom,prev_atom,MolID)

disp('Assuming all atom structs are unwrapped and whole...')

old_COM=zeros(numel(MolID),3);
for i=1:numel(MolID)
    temp = COM_atom(prev_atom([prev_atom.molid]==MolID(i)));
    position=COM;
    if size(new_atom,2) < 500
        temp_atom=COM_atom(new_atom); % This generates the COM position through assignin. You should use unwrapped molecules
        temp_atom=rmfield(temp_atom,'COM_x');
        temp_atom=rmfield(temp_atom,'COM_y');
        temp_atom=rmfield(temp_atom,'COM_z');
        temp_atom=rmfield(temp_atom,'element');
        temp_atom=rmfield(temp_atom,'Mw');
    else
        temp_atom = new_atom;
        COM=[mean([temp_atom.x]) mean([temp_atom.y]) mean([temp_atom.z])]; % Since COM_atom is a bit slow for large molecules, we do this for big molecules
    end

    temp_atom = translate_atom(new_atom,-COM+position,'all');
    [temp_atom.molid]=deal(MolID(i));

    ind_prev=find(ismember([prev_atom.molid],MolID(i)));
    % prev_atom(ind_prev)=[];
    % prev_atom = update_atom({prev_atom temp_atom});

    if MolID(i)==min([prev_atom.molid]) % Put the new molid first
        prev_atom=update_atom({temp_atom prev_atom(max(ind_prev)+1:end)});
    elseif MolID(i)==max([prev_atom.molid]) % Put the new molid last
        prev_atom=update_atom({temp_atom prev_atom(1:min(ind_prev)-1)});
    else % Put the new molid in between somewhere
        new_atom=update_atom({prev_atom(1:min(ind_prev)-1) temp_atom});
        prev_atom=update_atom({new_atom prev_atom(max(ind_prev)+1:end)});
    end

end

atom=prev_atom;

end

