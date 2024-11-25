%% order_attributes.m
% * This function order the struct attributes, or fields in a certain order.
%
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=order_attributes(atom) % Basic input arguments

function atom = order_attributes(atom)

if numel(fieldnames(atom))~=10
    defaultAttributes={'molid' 'molecule' 'resname' 'type' 'fftype' 'element' 'atnum' 'index' 'cn' 'neigh' 'bond' 'angle' 'x' 'y' 'z' 'vx' 'vy' 'vz' 'xfrac' 'yfrac' 'zfrac' 'mass' 'radius' 'Mw' 'COM_x' 'COM_y' 'COM_z' 'formalcharge' 'charge' 'bv' 'mean_bv' 'valence' 'Rdiff'  'occupancy' 'B'};
    atomAttributes=fieldnames(atom)';
    indDefault=find(ismember(defaultAttributes,atomAttributes));
    defaultAttributes=defaultAttributes(indDefault);
    ind_atom=find(ismember(atomAttributes,defaultAttributes));
    atomAttributes=atomAttributes(ind_atom);
    atom=orderfields(atom,unique({defaultAttributes{:} atomAttributes{:}},'stable'));
end


