%% translate_molid.m
% * This function translates the molid by a vector
%
%% Version
% 2.081
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = translate_molid(atom,[x y z],molid)

function atom = translate_molid(atom,trans_vec,molid)

disp('Translating MolID')

ind_MolID=ismember([atom.molid],molid);

x_shift=num2cell([atom(ind_MolID).x]+trans_vec(1)); [atom(ind_MolID).x]=deal(x_shift{:});

y_shift=num2cell([atom(ind_MolID).y]+trans_vec(2)); [atom(ind_MolID).y]=deal(y_shift{:});

z_shift=num2cell([atom(ind_MolID).z]+trans_vec(3)); [atom(ind_MolID).z]=deal(z_shift{:});

