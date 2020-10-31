%% merge_atom.m
% * This function returns the atom2 struct with atoms in the atom2 struct
% with a distance rmin [1x1 or 1x2] away from the atoms in the atom1 struct.
% * There is also a possibility to use a twin-range cutoff approach
% (suitable for OH2), by setting rmin(2) to a smaller value than rmin(1).
% * You must decide if atom1w should be wrapped or not before running this
% function
% * Do not wrap atomw1 here, do it before calling this function since
% Box1 does not always equal the full Box_dim, but rather a region in the
% full Box_dim
%
%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom2w = merge_atom(SOLUTE,Box_dim,SOLVENT,'index','C',1.4)
% # atom2w = merge_atom(SOLUTE,Box_dim,SOLVENT,'molid','Hw',[1.6 1.0])
%
function atom2w = merge_atom(atom1,Box1,atom2,type,Atom_label,rmin)
%%

% % % %% For testing the script
% % % clear all;
% % % format compact;
% % % atom1=import_atom('newTCH.gro');Box1=Box_dim;% atom1=translate_atom(atom1,[0 0 Box1(3)/2],'all');
% % % atom2=import_atom('spc216_1.00.gro');Box2=Box_dim;
% % % atom2 = replicate_atom(atom2,Box2,[5 5 3]);Box2=Box_dim;
% % % atom1=center_atom(atom1,Box1,'all','xyz');
% % % atom2=center_atom(atom2,Box2,'all','xyz');
% % % rmin=[2.2 3.2];
% % % limits=[0 0 0 Box1]
% % % Atom_label={'H' 'O'};
% % % type='molid';
% % %%

% if max([atom1.x])>Box1(1)
%     disp('Box1 smaller that Box_dim in x direction, perhaps using unwrapped atom1?')
%     disp('Will set Box1(1) to max([atom1.x])+.001')
%     Box1(1)=max([atom1.x])+.001;
% end
%
% if max([atom1.y])>Box1(2)
%     disp('Box1 smaller that Box_dim in y direction, perhaps using unwrapped atom1?')
%     disp('Will set Box1(2) to max([atom1.y])+.001')
%     Box1(2)=max([atom1.y])+.001;
% end
%
% if max([atom1.z])>Box1(3)
%     disp('Box1 smaller that Box_dim in z direction, perhaps using unwrapped atom1?')
%     disp('Will set Box1(3) to max([atom1.z])+.001')
%     Box1(3)=max([atom1.z])+.001;
% end

minAtom=Atom_label(1);
Box1max(1)=max([atom1.x]);
Box1max(2)=max([atom1.y]);
Box1max(3)=max([atom1.z]);
% Box1min(1)=min([atom1.x]);
% Box1min(2)=min([atom1.y]);
% Box1min(3)=min([atom1.z]);

if Box1max(1)>1.1*Box1(1) || Box1max(2)>1.1*Box1(2) || Box1max(3)>1.1*Box1(3)
    regional_solvation=1;
else
    regional_solvation=0;
end

atom1w=atom1; % Assuming atom1 is wrapped
% atom1w=wrap_atom(atom1,Box1);
% atom1w=atom1;% atom1 should hence be wrapped beforehand
atom2w=slice_molid(atom2,[0 0 0 Box1(1:3)]);
nAtoms1=size(atom1w,2);
nAtoms2=size(atom2w,2);
indvec=zeros(1,nAtoms2);

if size(atom2w,2)>6000
    % To do it stepwise to save memory...
    disp('Splitting up the distance matrix')
    size(atom1w);
    size(atom2w);
    X1_dist=pdist2([atom1w(1:floor(nAtoms1/2)).x]',[atom2w(1:floor(nAtoms2/2)).x]');
    Y1_dist=pdist2([atom1w(1:floor(nAtoms1/2)).y]',[atom2w(1:floor(nAtoms2/2)).y]');
    Z1_dist=pdist2([atom1w(1:floor(nAtoms1/2)).z]',[atom2w(1:floor(nAtoms2/2)).z]');
    X2_dist=pdist2([atom1w(1:floor(nAtoms1/2)).x]',[atom2w(floor((nAtoms2)/2)+1:end).x]');
    Y2_dist=pdist2([atom1w(1:floor(nAtoms1/2)).y]',[atom2w(floor((nAtoms2)/2)+1:end).y]');
    Z2_dist=pdist2([atom1w(1:floor(nAtoms1/2)).z]',[atom2w(floor((nAtoms2)/2)+1:end).z]');
    X3_dist=pdist2([atom1w(floor((nAtoms1)/2)+1:end).x]',[atom2w(1:floor(nAtoms2/2)).x]');
    Y3_dist=pdist2([atom1w(floor((nAtoms1)/2)+1:end).y]',[atom2w(1:floor(nAtoms2/2)).y]');
    Z3_dist=pdist2([atom1w(floor((nAtoms1)/2)+1:end).z]',[atom2w(1:floor(nAtoms2/2)).z]');
    X4_dist=pdist2([atom1w(floor((nAtoms1)/2)+1:end).x]',[atom2w(floor((nAtoms2)/2)+1:end).x]');
    Y4_dist=pdist2([atom1w(floor((nAtoms1)/2)+1:end).y]',[atom2w(floor((nAtoms2)/2)+1:end).y]');
    Z4_dist=pdist2([atom1w(floor((nAtoms1)/2)+1:end).z]',[atom2w(floor((nAtoms2)/2)+1:end).z]');
    
    if regional_solvation==0
        X1_dist(X1_dist>Box1(1)/2)=X1_dist(X1_dist>Box1(1)/2)-Box1(1);
        Y1_dist(Y1_dist>Box1(2)/2)=Y1_dist(Y1_dist>Box1(2)/2)-Box1(2);
        Z1_dist(Z1_dist>Box1(3)/2)=Z1_dist(Z1_dist>Box1(3)/2)-Box1(3);
        X2_dist(X2_dist>Box1(1)/2)=X2_dist(X2_dist>Box1(1)/2)-Box1(1);
        Y2_dist(Y2_dist>Box1(2)/2)=Y2_dist(Y2_dist>Box1(2)/2)-Box1(2);
        Z2_dist(Z2_dist>Box1(3)/2)=Z2_dist(Z2_dist>Box1(3)/2)-Box1(3);
        X3_dist(X3_dist>Box1(1)/2)=X3_dist(X3_dist>Box1(1)/2)-Box1(1);
        Y3_dist(Y3_dist>Box1(2)/2)=Y3_dist(Y3_dist>Box1(2)/2)-Box1(2);
        Z3_dist(Z3_dist>Box1(3)/2)=Z3_dist(Z3_dist>Box1(3)/2)-Box1(3);
        X4_dist(X4_dist>Box1(1)/2)=X4_dist(X4_dist>Box1(1)/2)-Box1(1);
        Y4_dist(Y4_dist>Box1(2)/2)=Y4_dist(Y4_dist>Box1(2)/2)-Box1(2);
        Z4_dist(Z4_dist>Box1(3)/2)=Z4_dist(Z4_dist>Box1(3)/2)-Box1(3);
    end
    
    X_dist=[X1_dist X2_dist;X3_dist X4_dist];
    Y_dist=[Y1_dist Y2_dist;Y3_dist Y4_dist];
    Z_dist=[Z1_dist Z2_dist;Z3_dist Z4_dist];
    
    dist_matrix=(X_dist.^2+Y_dist.^2+Z_dist.^2).^.5;
    clear X*dist Y*dist Z*dist;
elseif size(atom2w,2)>0
    % To do them all at once
    %disp('Single distance matrix')
    size([atom1w.x]');
    size([atom2w.x]');
    X_dist=pdist2([atom1w.x]',[atom2w.x]');
    Y_dist=pdist2([atom1w.y]',[atom2w.y]');
    Z_dist=pdist2([atom1w.z]',[atom2w.z]');
    
    if regional_solvation==0     
        X_dist(X_dist>Box1(1)/2)=X_dist(X_dist>Box1(1)/2)-Box1(1);
        Y_dist(Y_dist>Box1(2)/2)=Y_dist(Y_dist>Box1(2)/2)-Box1(2);
        Z_dist(Z_dist>Box1(3)/2)=Z_dist(Z_dist>Box1(3)/2)-Box1(3);
    end
    %%
    dist_matrix=(X_dist.^2+Y_dist.^2+Z_dist.^2).^.5;
    clear X*dist Y*dist Z*dist;
else
    disp('atom2w is empty')
end

% vmd([atom1w atom2w],Box1)
% disp('sizes X,Y,Z and distance matrix')
% size(X_dist)
%
% size(Y_dist)
%
% size(Z_dist)
%
% size(dist_matrix)



% This section loops though the atom1 atoms
% indvec=zeros(1,size(dist_matrix,2));
% for i=1:size(dist_matrix,1);
%     indvec(dist_matrix(i,:)<r)=1;
% end
if exist('dist_matrix','var')
    % This section loops though the atom2 atoms, using rsmall and rlarge, where rsmall operates on minAtom=Atom_label(1)
    if numel(rmin)==1
        rsmall=rmin;
        rlarge=rmin;
    else
        rsmall=rmin(1);
        rlarge=rmin(2);
    end
    
    for i=1:size(dist_matrix,2)
        if strncmpi([atom2w(i).type],minAtom,1)
            if any((dist_matrix(:,i)-rsmall)<0)
                indvec(1,i)=1;
            end
        else
            if any((dist_matrix(:,i)-rlarge)<0)
                indvec(1,i)=1;
            end
        end
    end
    
    clear dist_matrix;
    
    removed_molid=unique([atom2w(find(indvec)).molid]);
    removed_molid_index=[atom2w(ismember([atom2w.molid],removed_molid)).index];
    removed_index=find(indvec);
    
    pre_size=size(atom2w,2);
    if strcmpi(type,'molid')
        atom2w(ismember([atom2w.molid],removed_molid))=[];
    else
        atom2w(removed_index)=[];
    end
    post_size=size(atom2w,2);
    
    %     vmd([atom1w atom2w],Box1)
    
    % disp('Removed this many atoms')
    % pre_size-post_size
    
    assignin('caller','removed_molid',removed_molid);
    assignin('caller','removed_molid_index',removed_molid_index);
    assignin('caller','removed_index',removed_index);
    
end
% vmd([atom1w atom2w],Box1)




