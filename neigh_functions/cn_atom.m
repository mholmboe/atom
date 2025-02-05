%% cn_atom.m
% * This function tries to extract the coordination number of all the atom
% struct indexes and store it in the field [atom.cn]. The function calls
% either the bond_atom() or the cell_list_dist_matrix_atom() function
% directly to do this.
%
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=cn_atom(atom,Box_dim) % Basic input arguments
% # atom=cn_atom(atom,Box_dim,rmaxlong) % Allows setting the max cutoff

function atom = cn_atom(atom,Box_dim,varargin)

if nargin>2
    rmaxlong=varargin{1}; % Dummy value
else
    rmaxlong=2.25;
end

% if size(atom,2)>100000 && numel(Box_dim)<9
%     disp('Will use the cell list method')
%     disp('Does not work for triclinic systems...')
%     pause(5)
%     dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim,1.25,rmaxlong,rmaxlong,'more');
% else

% end

if nargin<4
    atom=bond_atom(atom,Box_dim,rmaxlong);
else
    Bond_index=varargin{2};
    atom=recalc_bond_atom(atom,Box_dim,Bond_index);
end

CN=num2cell(CoordNumber);
[atom.cn]=CN{:};

try
    assignin('caller','nBonds',nBonds);
    assignin('caller','radius_limit',radius_limit);
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','Neigh_index',Neigh_index);
    % assignin('caller','bond_matrix',dist_matrix);
    assignin('caller','dist_matrix',dist_matrix);
    assignin('caller','CoordNumber',CoordNumber);
    assignin('caller','Remove_ind',Remove_ind);
catch
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','CoordNumber',CoordNumber);
    assignin('caller','Remove_ind',Remove_ind);
    assignin('caller','dist_matrix',dist_matrix);
end

if nargin>3
    minCN=0;
    maxCN=max([atom.cn]);
    for i=minCN+1:maxCN+1
        assignin('caller',strcat('ind_CN',num2str(i-1)),find([atom.cn]==i-1));
    end
    Atom_labels=unique([atom.type]);
    for i=minCN+1:maxCN+1
        for a=1:length(Atom_labels)
            ind=find([atom.cn]==i-1);
            indtype=find(strcmp([atom.type],Atom_labels(a)));
            ind=intersect(ind,indtype);
            %             if ~isempty(ind)
            assignin('caller',strcat('ind_CN',num2str(i-1),'_',Atom_labels{a}),ind);
            %             end
        end
    end

end

atom = order_attributes(atom);

