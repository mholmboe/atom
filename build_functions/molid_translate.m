%% molid_translate.m
% * This function translates a certain MolID by the 1x3 translation vector
% * trans_vec. It's pretty straingh tforward.
%
%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = molid_translate(atom,[0 0 12],14)
%
function atom = molid_translate(atom,trans_vec,MolID)
%% 

disp('Translating MolID')

ind_MolID=ismember([atom.molid],MolID);

x_shift=num2cell([atom(ind_MolID).x]+trans_vec(1)); [atom(ind_MolID).x]=deal(x_shift{:});

y_shift=num2cell([atom(ind_MolID).y]+trans_vec(2)); [atom(ind_MolID).y]=deal(y_shift{:});

z_shift=num2cell([atom(ind_MolID).z]+trans_vec(3)); [atom(ind_MolID).z]=deal(z_shift{:});

