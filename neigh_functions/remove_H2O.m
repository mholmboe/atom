%% remove_H2O.m
% * This function removes H2O molecules, by searching for all atoms within
% rmin, which optionally can be set manually.
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=remove_H2O(atom,Box_dim) % Basic input arguments

function [atom,SOL] = remove_H2O(atom,Box_dim,varargin)

if nargin>2
    rmin=varargin{1};
else
    rmin=1.05;
end

dist_matrix = dist_matrix_atom(atom,Box_dim);

temp_atom=element_atom(atom);
O_ind=find(strncmp([temp_atom.type],{'O'},1));
H_ind=find(strncmp([temp_atom.type],{'H'},1));

SOL=[];
if numel(O_ind)>0
    rm_ind=[];n=1;molid=1;
    SOL=atom(O_ind(1));
    for i=1:length(O_ind)
        temp_ind=find(dist_matrix(O_ind(i),:)<rmin);
        if numel(temp_ind)>2
            rm_ind=[rm_ind temp_ind];
            ind=find(ismember(H_ind,temp_ind));
            if numel(ind)==2
                SOL(n)=atom(O_ind(i))
                SOL(n+1)=atom(H_ind(ind(1)));
                SOL(n+2)=atom(H_ind(ind(2)));
                [SOL(n:n+2).molid]=deal(molid);
                molid=molid+1;
                n=n+3;
            end
        end
    end
    numel(SOL)
    if size(SOL,2)>1
        [SOL.resname]=deal({'SOL'});
        SOL
        SOL=update_atom(SOL);
        SOL=bond_atom(SOL,Box_dim,rmin);
        SOL=update_atom(SOL);
        %     assignin('caller','SOL',SOL);
    else
        SOL=[];
    end
end

if size(SOL,2)<3
    SOL=[];
end

if numel(rm_ind)>0
    atom(unique(rm_ind))=[];
end
atom=update_atom(atom);

end

