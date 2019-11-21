%% wrap_molid.m
% * This function wraps the atom struct into the box
% * Todo, re-write according to the wrap_atom function...
%
%% Similar
% wrap_atom
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = wrap_molid(atom,Box_dim)

function atom = wrap_molid(atom,Box_dim)

atom = median_atom(atom);

ind_hix=find([atom.med_x]>=Box_dim(1));
x_shift=num2cell([[atom(ind_hix).x]-Box_dim(1)]');
[atom((ind_hix)).x]=deal(x_shift{:});

ind_hiy=find([atom.y]>=Box_dim(2));
y_shift=num2cell([[atom(ind_hiy).y]-Box_dim(2)]');
[atom((ind_hiy)).y]=deal(y_shift{:});

ind_hiz=find([atom.z]>=Box_dim(3));
z_shift=num2cell([[atom(ind_hiz).z]-Box_dim(3)]');
[atom((ind_hiz)).z]=deal(z_shift{:});

ind_lox=find([atom.med_x]<0);
x_shift=num2cell([[atom(ind_lox).x]+Box_dim(1)]');
[atom((ind_lox)).x]=deal(x_shift{:});

ind_loy=find([atom.med_y]<0);
y_shift=num2cell([[atom(ind_loy).y]+Box_dim(2)]');
[atom((ind_loy)).y]=deal(y_shift{:});

ind_loz=find([atom.med_z]<0);
z_shift=num2cell([[atom(ind_loz).z]+Box_dim(3)]');
[atom((ind_loz)).z]=deal(z_shift{:});

%assignin('caller','atom',atom);