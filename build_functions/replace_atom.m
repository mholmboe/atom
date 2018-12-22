%% replace_atom.m
% * This function replaces molid's in an atom struct with a new (single
% * molid) atom struct by placing the latters COM in the place of the former
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = replace_atom(in_atom,old_atom,molid_index)

function atom=replace_atom(new_atom,prev_atom,molid_index)

disp('Assuming all atom structs are unwrapped and whole...')

old_COM=zeros(numel(molid_index),3);
for i=1:numel(molid_index)
    temp = COM_atom(prev_atom([prev_atom.molid]==molid_index(i)));
    position=COM;
    if size(new_atom,2) < 500
        temp_atom=COM_atom(new_atom); % This generates the COM position through assignin. You should use unwrapped molecules
        temp_atom=rmfield(temp_atom,'COM_x');
        temp_atom=rmfield(temp_atom,'COM_y');
        temp_atom=rmfield(temp_atom,'COM_z');
        temp_atom=rmfield(temp_atom,'element');
        temp_atom=rmfield(temp_atom,'Mw');
    else
        COM=[mean([temp_atom.x]) mean([temp_atom.y]) mean([temp_atom.z])]; % Since COM_atom is a bit slow for large molecules, we do this for big molecules
    end
    
    temp_atom = translate_atom(new_atom,-COM+position,'all');
    [temp_atom.molid]=deal(molid_index(i));
    
    ind_prev=find(ismember([prev_atom.molid],molid_index(i)));
    if molid_index(i)==min([prev_atom.molid]) % Put the new molid first
        prev_atom=[temp_atom prev_atom(max(ind_prev)+1:end)];
    elseif molid_index(i)==max([prev_atom.molid]) % Put the new molid last
        prev_atom=[temp_atom prev_atom(1:min(ind_prev)-1)];
    else % Put the new molid in between somewhere
        prev_atom=[prev_atom(1:min(ind_prev)-1) temp_atom prev_atom(max(ind_prev)+1:end)];
    end

end

% Do we need this?
atom=update_atom(prev_atom);

