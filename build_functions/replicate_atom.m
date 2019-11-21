%% replicate_atom.m
% * This function replicates the atom struct and the orthogonal box dimensions
% * Triclinic version untestd but might work.
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = replicate_atom(atom,Box_dim,[6 4 1])
% # atom = replicate_atom(atom,Box_dim,[6 4 1],'yxz')
% * atom = replicate_atom(atom,Box_dim,[6 4 1],'xyz','addmolid')
%
function atom = replicate_atom(atom,Box_dim,replicate,varargin)

if numel(replicate)==1
    replicate(1)=replicate(1);
    replicate(2)=replicate(1);
    replicate(3)=replicate(1);
end
replicate(replicate==0)=1;

if numel(Box_dim)==1
    Box_dim(1)=Box_dim(1);
    Box_dim(2)=Box_dim(1);
    Box_dim(3)=Box_dim(1);
end

if length(Box_dim)==9
    Box_dim_tric=Box_dim;
    atom = orto_atom(atom,Box_dim_tric);
end

if nargin>3
    dim_order=char(varargin{1});
else
    dim_order='xyz';
end

if nargin>4
    addmolid=1;
else
    addmolid=0;
end

combinedatom_dim=atom;
molid=[atom(1).molid];
for j=1:3
    if j == strfind(dim_order,'x')
        if j == 1
            combinedatomx=atom;
            for i=1:replicate(1)
                if i > 1
                    newatom=atom;
                    if addmolid==1
                        molid=molid+1;
                        [newatom.molid]=deal(molid);
                    end
                    x_shift=num2cell([atom.x]+(i-1)*Box_dim(1)); [newatom.x]=deal(x_shift{:});
                    combinedatomx=[combinedatomx newatom];
                end
            end
        else
            combinedatomx=combinedatom_dim;
            for i=1:replicate(1)
                if i > 1
                    newatom=combinedatom_dim;
                    x_shift=num2cell([combinedatom_dim.x]+(i-1)*Box_dim(1)); [newatom.x]=deal(x_shift{:});
                    combinedatomx=[combinedatomx newatom];
                end
            end
        end
        combinedatom_dim=combinedatomx;
        disp('x - dim')
        length(combinedatom_dim)
        
    elseif j == strfind(dim_order,'y')
        if j == 1
            combinedatomy=atom;
            for i=1:replicate(2)
                if i > 1
                    newatom=atom;
                    if addmolid==1
                        molid=molid+1;
                        [newatom.molid]=deal(molid);
                    end
                    y_shift=num2cell([atom.y]+(i-1)*Box_dim(2)); [newatom.y]=deal(y_shift{:});
                    combinedatomy=[combinedatomy newatom];
                end
            end
        else
            combinedatomy=combinedatom_dim;
            for i=1:replicate(2)
                if i > 1
                    newatom=combinedatom_dim;
                    y_shift=num2cell([combinedatom_dim.y]+(i-1)*Box_dim(2)); [newatom.y]=deal(y_shift{:});
                    combinedatomy=[combinedatomy newatom];
                end
            end
        end
        combinedatom_dim=combinedatomy;
        disp('y - dim')
        length(combinedatom_dim)
        
    elseif j == strfind(dim_order,'z')
        if j == 1
            combinedatomz=atom;
            for i=1:replicate(3)
                if i > 1
                    newatom=atom;
                    if addmolid==1
                        molid=molid+1;
                        [newatom.molid]=deal(molid);
                    end
                    z_shift=num2cell([atom.z]+(i-1)*Box_dim(3)); [newatom.z]=deal(z_shift{:});
                    combinedatomz=[combinedatomz newatom];
                end
            end
        else
            combinedatomz=combinedatom_dim;
            for i=1:replicate(3)
                if i > 1
                    newatom=combinedatom_dim;
                    z_shift=num2cell([combinedatom_dim.z]+(i-1)*Box_dim(3)); [newatom.z]=deal(z_shift{:});
                    combinedatomz=[combinedatomz newatom];
                end
            end
        end
        combinedatom_dim=combinedatomz;
        disp('z - dim')
        length(combinedatom_dim)
    end
end

atom=combinedatom_dim;

Box_dim=[Box_dim(1)*replicate(1) Box_dim(2)*replicate(2) Box_dim(3)*replicate(3)];

if exist('Box_dim_tric','var')
    atom = triclinic_atom(atom,Box_dim,[replicate(1)*Box_dim_tric(6) replicate(3)*Box_dim_tric(8) replicate(2)*Box_dim_tric(9)],'tiltfactors');
    Box_dim=triclinic_Box_dim;
end

% combinedatomx=atom;
% for i=1:replicate(1);
%     if i > 1;
%         newatom=atom;
%         x_shift=num2cell([atom.x]+(i-1)*Box_dim(1)); [newatom.x]=deal(x_shift{:});
%         combinedatomx=[combinedatomx newatom];
%     end
% end
% length(combinedatomx);
%
% combinedatomy=combinedatomx;
% for i=1:replicate(2);
%     if i > 1;
%         newatom=combinedatomx;
%         y_shift=num2cell([combinedatomx.y]+(i-1)*Box_dim(2)); [newatom.y]=deal(y_shift{:});
%         combinedatomy=[combinedatomy newatom];
%     end
% end
% length(combinedatomy);
%
% combinedatomz=combinedatomy;
% for i=1:replicate(3);
%     if i > 1;
%         newatom=combinedatomy;
%         z_shift=num2cell([combinedatomy.z]+(i-1)*Box_dim(3)); [newatom.z]=deal(z_shift{:});
%         combinedatomz=[combinedatomz newatom];
%     end
% end
% length(combinedatomz);
%
% Box_dim=[Box_dim(1)*replicate(1) Box_dim(2)*replicate(2) Box_dim(3)*replicate(3)];
%
% atom=combinedatomz;

atom=update_atom(atom);

XYZ_data=[[atom.x]' [atom.y]' [atom.z]']; XYZ_labels=[atom.type]';

assignin('caller','Box_dim',Box_dim);
assignin('caller','XYZ_labels',XYZ_labels);
assignin('caller','XYZ_data',XYZ_data);