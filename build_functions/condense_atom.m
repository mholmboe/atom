%% condense_atom.m
% * This function tries to minimize the box size and remove gaps between molecules along x,y,z
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = wrap_atom_func(atom,Box_dim);
%
function atom = condense_atom(atom,Box_dim,s)

x_shift=num2cell([[atom.x]-min([atom.x])]'); [atom(:).x]=deal(x_shift{:});
y_shift=num2cell([[atom.y]-min([atom.y])]'); [atom(:).y]=deal(y_shift{:});
z_shift=num2cell([[atom.z]-min([atom.z])]'); [atom(:).z]=deal(z_shift{:});

Box_dim(1)=max([atom.x]);Box_dim(2)=max([atom.y]);Box_dim(3)=max([atom.z]);

for repeat=1:10
    %s=0.5;
    x_shift=0;y_shift=0;z_shift=0;
    for i=1:s:Box_dim(1)
        j=i-s;
        ind_hi=find([atom.x]>=(j));
        ind_lo=find([atom.x]<=(i));
        ind=intersect(ind_lo,ind_hi);
        if length(ind)==0 && length(ind_hi) > 0
            shift=num2cell([[atom(ind_hi).x]-s]');
            [atom(ind_hi).x]=deal(shift{:});
            x_shift=x_shift+s;
        end
    end
    Box_dim(1)=Box_dim(1)-x_shift
    for i=1:s:Box_dim(2)
        j=i-s;
        ind_hi=find([atom.y]>=(j));
        ind_lo=find([atom.y]<=(i));
        ind=intersect(ind_lo,ind_hi);
        if length(ind)==0 && length(ind_hi) > 0
            shift=num2cell([[atom(ind_hi).y]-s]');
            [atom(ind_hi).y]=deal(shift{:});
            y_shift=y_shift+s;
        end
    end
    Box_dim(2)=Box_dim(2)-y_shift
    for i=1:s:Box_dim(3)
        j=i-s;
        ind_hi=find([atom.z]>=(j));
        ind_lo=find([atom.z]<=(i));
        ind=intersect(ind_lo,ind_hi);
        if length(ind)==0 && length(ind_hi) > 0
            shift=num2cell([[atom(ind_hi).z]-s]');
            [atom(ind_hi).z]=deal(shift{:});
            z_shift=z_shift+s;
        end
    end
    Box_dim(3)=Box_dim(3)-z_shift
end

assignin('caller','Box_dim',Box_dim)
