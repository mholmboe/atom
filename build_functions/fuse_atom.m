%% fuse_atom.m
% * This function tries to fuse all sites within a certain radii, rmax.
%
%% Version
% 2.09
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = fuse_atom(atom,Box_dim) % Basic input arguments, rmax set to 0.85
% * atom = fuse_atom(atom,Box_dim,1.0) % rmax set to 1 Å

function atom = fuse_atom(atom,Box_dim,varargin) % ,rmax);

fused_atom=atom;

dist_matrix=dist_matrix_atom(atom,Box_dim);

if nargin==2
    rmax=0.85;
else
    rmax=varargin{1};
end

i=1;rmind_tot=[];
while i < size(fused_atom,2)
    rmind=find(dist_matrix(:,i)<rmax)';
    if numel(rmind)>1
        x1=[fused_atom(i).x];
        y1=[fused_atom(i).y];
        z1=[fused_atom(i).z];
        
%         fused_atom(rmind) = translate_atom(fused_atom(rmind),[Box_dim(1)/2-x1 Box_dim(2)/2-y1 Box_dim(3)/2-z1]);
%         fused_atom(rmind) = wrap_atom(fused_atom(rmind),Box_dim);

        [fused_atom(i).x]=fused_atom(i).x-mean(X_dist(rmind,i));% mean([fused_atom(rmind).x]);
        [fused_atom(i).y]=fused_atom(i).y-mean(Y_dist(rmind,i));% mean([fused_atom(rmind).y]);
        [fused_atom(i).z]=fused_atom(i).z-mean(Z_dist(rmind,i));% mean([fused_atom(rmind).z]);
        
%         fused_atom(rmind) = translate_atom(fused_atom(rmind),[-Box_dim(1)/2+x1 -Box_dim(2)/2+y1 -Box_dim(3)/2+z1]);
        rmind_tot=[rmind_tot rmind(rmind>i)];
    end
    i=i+1;
    
    if mod(i,100)==1
        if i > 1
            i-1
        end
    end
    
end
fused_atom(rmind_tot)=[];

atom=update_atom(fused_atom);

% if isstruct(fused_atom)
%     try
%         if ~isfield(atom,'element')
%             fused_atom=rmfield(fused_atom,'element');
%         end
%     catch
%     end
%
% %     try
% %         if isfield(fused_atom,'xfrac')
% %             fused_atom=rmfield(fused_atom,'xfrac');
% %             fused_atom=rmfield(fused_atom,'yfrac');
% %             fused_atom=rmfield(fused_atom,'zfrac');
% %         end
% %     catch
% %     end
% end
