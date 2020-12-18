%% radius_atom.m
% * This function fetches the ion radius from clayff or interface or interface2015 ff's and
% adds a new radius field to the atoms struct
%
%% Version
% 2.081
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = radius_atom(atom,clayff,'spc')
% # atom = radius_atom(atom,interface,'tip3p')

function atom = radius_atom(atom,ffname,watermodel)
%% 

if strcmpi(ffname,'clayff')
    clayff_param(unique([atom.type]),watermodel);
    % Check the charge after AssignClayff.m
    for i=1:length(atom)
        if strncmpi([atom(i).type],{'Hw'},2)
            ind=strncmpi({'Hw'},[forcefield.clayff.type],2);
        else
            ind=strcmpi([atom(i).type],[forcefield.clayff.type]);
        end
        atom(i).radius=[forcefield.clayff(ind).radius];
    end
%     atom = charge_interface_atom(atom,Box_dim,varargin);
elseif strcmpi(ffname,'interface')
    interface_param(unique([atom.type]),watermodel);
    for i=1:length(atom)
        if strncmpi([atom(i).type],{'Hw'},2)
            ind=strncmpi({'Hw'},[forcefield.interface.type],2);
        else
            ind=strcmp([atom(i).type],[forcefield.interface.type]);
        end
        atom(i).radius=[forcefield.interface(ind).radius];
    end
elseif strcmpi(ffname,'interface15')
    interface15_param(unique([atom.type]),watermodel);
    for i=1:length(atom)
        if strncmpi([atom(i).type],{'Hw'},2)
            ind=strncmpi({'Hw'},[forcefield.interface15.type],2);
        else
            ind=strcmp([atom(i).type],[forcefield.interface15.type]);
        end
        atom(i).radius=[forcefield.interface15(ind).radius];
    end
end
