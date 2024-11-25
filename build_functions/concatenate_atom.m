%% concatenate_atom.m
% * This old function concatenats atom sections. Use update_atom({atom_1 atom_2}) instead...
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = concatenate_atom(atom_1,atom_2) % Basic input arguments
%
function atom = concatenate_atom(atom,varargin)
disp('This old function concatenats atom sections. Use update_atom({atom_1 atom_2}) instead...')

% atom=[atom_1 atom_2]; % This also works in principle, but does not update any indexes

if nargin>1
    for i=1:nargin-1
        atom=update_atom({atom varargin{i}});
    end
end

end