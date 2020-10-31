%% position_molid.m
% * This function movies a molid (COM) to a certain position
%
%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = position_molid(atom,[0 0 8],1)

function atom = position_molid(atom,position_vec,MolID)

disp('Moving MolID (COM) to certain position')

% atom = median_atom_func(atom);

atom = COM_atom(atom,MolID);

ind=ismember([atom.molid],MolID);

x_shift=num2cell([atom(ind).x]-[atom(ind).COM_x]+position_vec(1)); [atom(ind).x]=deal(x_shift{:});

y_shift=num2cell([atom(ind).y]-[atom(ind).COM_y]+position_vec(2)); [atom(ind).y]=deal(y_shift{:});

z_shift=num2cell([atom(ind).z]-[atom(ind).COM_z]+position_vec(3)); [atom(ind).z]=deal(z_shift{:});
