%% unwrap_atom.m
% * This function unwraps the atom struct along the dimension dim
% * Tested 21/07/2016, there has been bugs.. does it work?
%
%% Version
% 2.06
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = unwrap_atom(atom,Box_dim,'x')
% # atom = unwrap_atom(atom,Box_dim,'xyz')
%
function atom = unwrap_atom(atom,Box_dim,dim)

atom_Tot=atom;

for m=unique([atom_Tot.molid])

atom = atom_Tot([atom_Tot.molid]==m);
%Broken_molid=[];
if sum(ismember(dim,'x'))>0
    Broken_molid=[];
    disp('Unwrap in x')
    
%     for i=1:max([atom(:).molid])
%         ind=find([atom.molid]==i);
%         for j=ind(2:end)
%             if abs(atom(j).x-atom(j-1).x)>Box_dim(1)/2;
%                 Broken_molid=[Broken_molid i];
%             end
%         end
%     end
%     
%     Broken_ind=find(ismember([atom.molid],Broken_molid));
%     
%     atom = median_atom(atom);
%     
%     ind_med_hi=find([atom.med_x]>=Box_dim(1)/2);
%     ind_med_hi=intersect(ind_med_hi,Broken_ind);
%     ind_lo=find([atom.x]<Box_dim(1)/2&[atom.x]<([atom.med_x]-1.1*Box_dim(1)/2));
%     ind_wrap=intersect(ind_med_hi,ind_lo);
%     shift=num2cell([atom(ind_wrap).x]+Box_dim(1));
%     [atom((ind_wrap)).x]=deal(shift{:});
%     
%     ind_med_lo=find([atom.med_x]<=Box_dim(1)/2);
%     ind_med_lo=intersect(ind_med_lo,Broken_ind);
%     ind_hi=find([atom.x]>Box_dim(1)/2);
%     ind_hi=find([atom.x]>Box_dim(1)/2&[atom.x]>([atom.med_x]+1.1*Box_dim(1)/2));
%     ind_wrap=intersect(ind_med_lo,ind_hi);
%     shift=num2cell([atom(ind_wrap).x]-Box_dim(1));
%     [atom((ind_wrap)).x]=deal(shift{:});
%     
%     ind_2hi=find([atom.x]>Box_dim(1)/2+.7*Box_dim(1));
%     shift=num2cell([atom(ind_2hi).x]-Box_dim(1));
%     [atom((ind_2hi)).x]=deal(shift{:});
%     
%     ind_2lo=find([atom.x]<Box_dim(1)/2-.7*Box_dim(1));
%     shift=num2cell([atom(ind_2lo).x]+Box_dim(1));
%     [atom((ind_2lo)).x]=deal(shift{:});
    
        for k=2:size(atom,2)
            if (atom(k).x-mean([atom(1:k-1).x]))<-Box_dim(1)/2
                atom(k).x=atom(k).x+Box_dim(1);
            elseif (atom(k).x-mean([atom(1:k-1).x]))>Box_dim(1)/2
                atom(k).x=atom(k).x-Box_dim(1);
            end
        end
    
end

if sum(ismember(dim,'y'))>0
    Broken_molid=[];
    disp('Unwrap in y')
    
%     for i=1:max([atom(:).molid])
%         ind=find([atom.molid]==i);
%         for j=ind(2:end)
%             if abs(atom(j).y-atom(j-1).y)>Box_dim(2)/2;
%                 Broken_molid=[Broken_molid i];
%             end
%         end
%     end
%     
%     Broken_ind=find(ismember([atom.molid],Broken_molid));
%     
%     ind_med_hi=find([atom.med_y]>=Box_dim(2)/2);
%     ind_med_hi=intersect(ind_med_hi,Broken_ind);
%     ind_lo=find([atom.y]<Box_dim(2)/2&[atom.y]<([atom.med_y]-1.1*Box_dim(2)/2));
%     ind_wrap=intersect(ind_med_hi,ind_lo);
%     shift=num2cell([atom(ind_wrap).y]+Box_dim(2));
%     [atom((ind_wrap)).y]=deal(shift{:});
%     
%     ind_med_lo=find([atom.med_y]<=Box_dim(2)/2);
%     ind_med_lo=intersect(ind_med_lo,Broken_ind);
%     ind_hi=find([atom.y]>Box_dim(2)/2&[atom.y]>([atom.med_y]+1.1*Box_dim(2)/2));
%     ind_wrap=intersect(ind_med_lo,ind_hi);
%     shift=num2cell([atom(ind_wrap).y]-Box_dim(2));
%     [atom((ind_wrap)).y]=deal(shift{:});
%     
%     ind_2hi=find([atom.y]>Box_dim(2)/2+.7*Box_dim(2))
%     shift=num2cell([atom(ind_2hi).y]-Box_dim(2));
%     [atom((ind_2hi)).y]=deal(shift{:});
%     
%     ind_2lo=find([atom.y]<Box_dim(2)/2-.7*Box_dim(2))
%     shift=num2cell([atom(ind_2lo).y]+Box_dim(2));
%     [atom((ind_2lo)).y]=deal(shift{:});
    
        for k=2:size(atom,2)
            if (atom(k).y-mean([atom(1:k-1).y]))<-Box_dim(2)/2
                atom(k).y=atom(k).y+Box_dim(2);
            elseif (atom(k).y-mean([atom(1:k-1).y]))>Box_dim(2)/2
                atom(k).y=atom(k).y-Box_dim(2);
            end
        end
    
end

if sum(ismember(dim,'z'))>0
    Broken_molid=[];
    disp('Unwrap in z')
    
%     for i=1:max([atom(:).molid])
%         ind=find([atom.molid]==i);
%         for j=ind(2:end)
%             if abs(atom(j).z-atom(j-1).z)>Box_dim(3)/2;
%                 Broken_molid=[Broken_molid i];
%             end
%         end
%     end
%     
%     Broken_ind=find(ismember([atom.molid],Broken_molid));
%     
%     ind_med_hi=find([atom.med_z]>=Box_dim(3)/2);
%     ind_med_hi=intersect(ind_med_hi,Broken_ind);
%     ind_lo=find([atom.z]<Box_dim(3)/2&[atom.z]<([atom.med_z]-1.1*Box_dim(3)/2));
%     ind_wrap=intersect(ind_med_hi,ind_lo);
%     shift=num2cell([atom(ind_wrap).z]+Box_dim(3));
%     [atom((ind_wrap)).z]=deal(shift{:});
%     
%     ind_med_lo=find([atom.med_z]<=Box_dim(3)/2);
%     ind_med_lo=intersect(ind_med_lo,Broken_ind);
%     ind_hi=find([atom.z]>Box_dim(3)/2&[atom.z]>([atom.med_z]+1.1*Box_dim(3)/2));
%     ind_wrap=intersect(ind_med_lo,ind_hi);
%     shift=num2cell([atom(ind_wrap).z]-Box_dim(3));
%     [atom((ind_wrap)).z]=deal(shift{:});
%     
%     ind_2hi=find([atom.z]>Box_dim(3)/2+.7*Box_dim(3));
%     shift=num2cell([atom(ind_2hi).z]-Box_dim(3));
%     [atom((ind_2hi)).z]=deal(shift{:});
%     
%     ind_2lo=find([atom.z]<Box_dim(3)/2-.7*Box_dim(3));
%     shift=num2cell([atom(ind_2lo).z]+Box_dim(3));
%     [atom((ind_2lo)).z]=deal(shift{:});

        for k=2:size(atom,2)
            if (atom(k).z-mean([atom(1:k-1).z]))<-Box_dim(3)/2
                atom(k).z=atom(k).z+Box_dim(3);
            elseif (atom(k).z-mean([atom(1:k-1).z]))>Box_dim(3)/2
                atom(k).z=atom(k).z-Box_dim(3);
            end
        end
        
        atom = median_atom(atom);
        
        xshift=0;yshift=0;zshift=0;
        if [atom(1).med_x]>Box_dim(1)
            xshift=-Box_dim(1);
        elseif [atom(1).med_x]<0
            xshift=Box_dim(1);
        end
        
        if [atom(1).med_y]>Box_dim(2)
            yshift=-Box_dim(2);
        elseif [atom(1).med_y]<0
            yshift=Box_dim(2);
        end
        
        if [atom(1).med_z]>Box_dim(3)
            zshift=-Box_dim(3);
        elseif [atom(1).med_z]<0
            zshift=Box_dim(3);
        end
    
        atom=translate_atom(atom,[xshift yshift zshift],'all');
        
        
end

[atom_Tot([atom_Tot.molid]==m).x]=atom.x;
[atom_Tot([atom_Tot.molid]==m).y]=atom.y;
[atom_Tot([atom_Tot.molid]==m).z]=atom.z;

end

% atom_Tot = rmfield(atom_Tot,'med_x');
% atom_Tot = rmfield(atom_Tot,'med_y');
% atom_Tot = rmfield(atom_Tot,'med_z');

atom = atom_Tot;




