clear all;
format short;
%% This scripts allows us to take 1 unit cell of something and perform isomorphous substitution

%% Edit this section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UCinX,UCinY,UCinY will not be used when building edges
UCinX = 6;                      % Number of Unit cells in x direction
UCinY = 4;                      % Number of Unit cells in y direction
UCinZ = 1;                      % Number of Unit cells in z direction. This must be 1 (I think)
NumOctSubst=16;                 % How many octahedral substituitions do you want,
NumTetSubst=0;                  % How many tetrahedral substituitions do you want
O1={'Al'}; O2={'Mgo'};  % Mgo   % O1 will get replaced by O2 NumOctSubst times
T1={'Si'}; T2={'Alt'};  % Alt   % T1 will get replaced by T2 NumTetSubst times
O2O2_dist=4.5;                  % Minimum O2/O2 substitutions distance in Å, you may decrease it from 4.5 if you add alot of charge
T2T2_dist=4.5;                  % Minimum T2/T2 substitutions distance in Å, you may decrease it from 4.5 if you add alot of charge

filename='1xMMT.gro';        % Coordinates taken from Bickmore 2003 and made non-triclininc
                                              % UC dimension is 0.51488   0.899790   0.98409 in nm! all angles are 90deg
filename_out='MMT_4.gro';  % Set the name of the output file

%% Do not edit anything below unless you know what you do...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import the input gro structure file
import_atom(filename); disp('Found .gro file');

% Replicate the unit cell
atom = replicate_atom(atom,Box_dim,[UCinX UCinY UCinZ]);
atom = center_atom(atom,[Box_dim(1) Box_dim(2) 2*Box_dim(3)],'all','z');
atom = translate_atom(atom,[0 0 -median([atom.z])],'all');

% Limits for the substitutions, only needed when building edges
lolimit=-1000; % Arbitrary low number
hilimit= 1000;  % Arbitrary high number
dimension=2;   %x=1, y=2, z=3

% Total n of substitutions
NumTotalSubst=NumOctSubst+NumTetSubst;

XYZ_labels=[atom.type]';
XYZ_data = [[atom.x]' [atom.y]' [atom.z]'];

if NumOctSubst>0;
    
    atom=bond_angle_atom(atom,Box_dim,O2O2_dist,O2O2_dist);
    O1_Index=find(strcmpi(O1,strtrim(XYZ_labels(:,1))));
    O1_labels=XYZ_labels(O1_Index);
    O1_data=XYZ_data(O1_Index,1:3);
    Ave_Oct_z=mean(O1_data(:,3));
    rand_O1_Index=O1_Index(randperm(length(O1_Index)));

    i=2; nOctlo=0; nOcthi=0; nOctmid=0; Oct_subst_index=[];%rand_O1_Index(1);
    while (nOctlo+nOcthi+nOctmid)<NumOctSubst;
        [A1,A2]=find((ismember(Angle_index,Oct_subst_index)));
        [B1,B2]=find((ismember(Angle_index,rand_O1_Index(i))));
        C=intersect(A1,B1,'rows');
        [D1,D2]=find((ismember(Bond_index,Oct_subst_index)));
        [E1,E2]=find((ismember(Bond_index,rand_O1_Index(i))));
        F=intersect(D1,E1,'rows');
        if length(C)<1 && length(F)<1 && nOctlo < NumOctSubst && lolimit < XYZ_data(rand_O1_Index(i),dimension) && hilimit > XYZ_data(rand_O1_Index(i),dimension);
            if nOctlo < NumOctSubst/2 && XYZ_data(rand_O1_Index(i),3)<Ave_Oct_z; %&& ceil(XYZ_data(rand_O1_Index(i),3)*100)/100<=Ave_Oct_z;
                Oct_subst_index=[Oct_subst_index; rand_O1_Index(i)];
                nOctlo=nOctlo+1;
            elseif nOcthi < NumOctSubst/2 && XYZ_data(rand_O1_Index(i),3)>Ave_Oct_z; %&& ceil(XYZ_data(rand_O1_Index(i),3)*100)/100>=Ave_Oct_z;
                Oct_subst_index=[Oct_subst_index; rand_O1_Index(i)];
                nOcthi=nOcthi+1;
            elseif (nOctlo+nOcthi+nOctmid) < NumOctSubst && XYZ_data(rand_O1_Index(i),3)==Ave_Oct_z;
                Oct_subst_index=[Oct_subst_index; rand_O1_Index(i)];
                nOctmid=nOctmid+1;
            end
            
        end
        
        if i == length(O1_data);
            break
            disp('Stopped the loop')
        end
        i=i+1;
    end
    XYZ_labels(Oct_subst_index)=O2;
    O2_labels=XYZ_labels(Oct_subst_index);
    O2_data=XYZ_data(Oct_subst_index,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NumTetSubst>0;
    
    atom=bond_angle_atom(atom,Box_dim,T2T2_dist,T2T2_dist);
    T1_Index=find(strcmpi(T1,strtrim(XYZ_labels(:,1))));
    T1_labels=XYZ_labels(T1_Index);
    T1_data=XYZ_data(T1_Index,1:3);
    Ave_Tet_z=mean(T1_data(:,3));
    rand_T1_Index=T1_Index(randperm(length(T1_Index)));

    i=2; nTetlo=0; nTethi=0; All_subst_index=Oct_subst_index; Tet_subst_index=[];%(1)= rand_T1_Index(1);
    while (nTetlo+nTethi)<NumTetSubst;
        [A1,A2]=find((ismember(Angle_index,All_subst_index)));
        [B1,B2]=find((ismember(Angle_index,rand_T1_Index(i))));
        C=intersect(A1,B1,'rows');
        [D1,D2]=find((ismember(Bond_index,All_subst_index)));
        [E1,E2]=find((ismember(Bond_index,rand_T1_Index(i))));
        F=intersect(D1,E1,'rows');
        if length(C)<1 && length(F)<1 && nTetlo < NumTetSubst && lolimit < XYZ_data(rand_T1_Index(i),dimension) && hilimit > XYZ_data(rand_T1_Index(i),dimension);
            if nTetlo < NumTetSubst/2 && ceil(XYZ_data(rand_T1_Index(i),3)*100)/100<=Ave_Tet_z ;
                All_subst_index=[All_subst_index; rand_T1_Index(i)];
                Tet_subst_index=[Tet_subst_index; rand_T1_Index(i)];
                nTetlo=nTetlo+1;
            elseif nTethi < NumTetSubst/2 && ceil(XYZ_data(rand_T1_Index(i),3)*100)/100>=Ave_Tet_z;
                All_subst_index=[All_subst_index; rand_T1_Index(i)];
                Tet_subst_index=[Tet_subst_index; rand_T1_Index(i)];
                nTethi=nTethi+1;
            end
        end
        if i == length(T1_data);
            break
            disp('Stopped the loop')
        end
        i=i+1;
    end
    
    XYZ_labels(Tet_subst_index)=T2;
    T2_labels=XYZ_labels(Tet_subst_index);
    T2_data=XYZ_data(Tet_subst_index,:);
    
    if nTetlo==nTethi && (nTetlo+nTethi) == NumTetSubst;
        disp('Tetrahedral substitution success!!!')
    end
    
end

if nOctlo==nOcthi && (nOctlo+nOcthi) == NumOctSubst;
    disp('Octahedral substitution success!!!')
end

nAtoms=size(XYZ_data,1);
for i=1:nAtoms;
    atom(i).type   = XYZ_labels(i,1);
    atom(i).fftype = XYZ_labels(i,1);
    atom(i).x      = XYZ_data(i,1);
    atom(i).y      = XYZ_data(i,2);
    atom(i).z      = XYZ_data(i,3);
end

if NumOctSubst>0;
atom_O2=atom(find(strcmpi([atom.type],O2)));
O2_distmatrix=dist_matrix_atom(atom_O2,Box_dim);
disp('Minimum O2O2_dist is in Å')
min(O2_distmatrix(2:end,1))
end

if NumTetSubst>0;
atom_T2=atom(find(strcmpi([atom.type],T2)));
T2_distmatrix=dist_matrix_atom(atom_T2,Box_dim);
disp('Minimum T2T2_dist is in Å')
min(T2_distmatrix(2:end,1))
end

%%%%%% Write structure to file %%%%%%
write_atom_gro(atom,Box_dim,filename_out);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Plot the new structure in VMD %%%%%%
vmd(atom(ismember([atom.type],[O1 O2 T2])),Box_dim)
%%%%%%%%%
