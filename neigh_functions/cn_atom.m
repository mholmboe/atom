%% cn_atom.m
% * This function tries to extract the coordination number of all the atom
% struct indexes and store it in the field atom.cn. The function calls
% either the bond_atom() or the cell_list_dist_matrix_atom() function
% directly to do this.
%
%
%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=cn_atom(atom,Box_dim)
% # atom=cn_atom(atom,Box_dim,2.1) % Last argument is the max bond distance

function atom = cn_atom(atom,Box_dim,varargin)

if nargin>2
    rmaxlong=varargin{1}; % Dummy value
else
    rmaxlong=2.25;
end

if size(atom,2)<15000
    atom=bond_atom(atom,Box_dim,rmaxlong);
else
    dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim,1.25,rmaxlong,rmaxlong,'more');
end

CN=num2cell(CoordNumber);
[atom.cn]=CN{:};


assignin('caller','Bond_index',Bond_index);
assignin('caller','CoordNumber',CoordNumber);
assignin('caller','Remove_ind',Remove_ind);

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



