%% substitute_atom.m
% * This scripts performs isomorphous substitution, by replacing some O1->O2 atomtypes and optionally T1->T2 atomtypes
% * varargin should be something like ,NumTetrahedralSubst,T1,T2,minT2T2_dist)
% * aotm is the atom struct
% * Box_dimensions is the 
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom = substitute_atom(atom,Box_dim,5,'Ca','Mgo',5.5)
% * atom = substitute_atom(atom,Box_dim,5,'Ca','Mgo',5.5,2,'Si','Al',5.5)


function atom = substitute_atom(atom,Box_dim,NumOctSubst,O1,O2,minO2O2_dist,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you want to replicate the box, edit this section. UCinX,UCinY,UCinY will not be used when building edges
% clear all;
% UCinX = 6;                      % Number of Unit cells in x direction
% UCinY = 4;                      % Number of Unit cells in y direction
% UCinZ = 1;                      % Number of Unit cells in z direction. This must be 1 (I think)
% % %% Some suggested settings
% NumOctSubst=16;                 % How many octahedral substituitions do you want,
% NumTetSubst=4;                  % How many tetrahedral substituitions do you want
% O1={'Al'}; O2={'Mgo'};  % Mgo   % O1 will get replaced by O2 NumOctSubst times
% T1={'Si'}; T2={'Alt'};  % Alt   % T1 will get replaced by T2 NumTetSubst times
% minO2O2_dist=5.5;                  % Minimum O2/O2 substitutions distance in Å, you may decrease it to 4.5 if you add alot of charge
% minT2T2_dist=5.5;                  % Minimum T2/T2 substitutions distance in Å, you may decrease it to 4.5 if you add alot of charge


if ~iscell(O1);
    disp('Converting O1 to cell')
    O1={O1};
end
if ~iscell(O2);
    disp('Converting O2 to cell')
    O2={O2};
end

if nargin > 6;
    NumTetSubst=cell2mat(varargin(1));
    if iscell(varargin(2))>0
%         T1=varargin{2}(:)
%         T2=varargin{3}(:)
        T1=varargin{2}
        T2=varargin{3}
    else
        T1=varargin(2)
        T2=varargin(3)
    end
    minT2T2_dist=cell2mat(varargin(4))
    
    if ~iscell(T1);
        disp('Converting T1 to cell')
        T1={T1}
    end
    
    if ~iscell(T2);
        disp('Converting T2 to cell')
        T2={T2}
    end
else
    NumTetSubst=0;
    T1=O1; T2=O2;  % Alt   % T1 will get replaced by T2 NumTetSubst times
    minT2T2_dist=5.5;      % Minimum T2/T2 substitutions distance in Å, you may decrease it to 4.5 if you add alot of charge
end

% filename='1xpyro_skipper.gro';        % Coordinates taken from Bickmore 2003 and made non-triclininc, UC dimension is 0.51488   0.899790   0.98409 in nm! all angles are 90deg
% filename_out='MMT_6_interface.gro';  % Set the name of the output file
% filenameout='MMT1';
% atom=import_atom('1xMMT.gro');
% %% Do not edit anything below unless you know what you do...
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Replicate the unit cell
% atom = replicate_atom(atom,Box_dim,[6 4 1]);
% atom = center_atom(atom,[Box_dim(1) Box_dim(2) 2*Box_dim(3)],'all','z');
% atom = translate_atom(atom,[0 0 -median([atom.z])],'all');

% Limits for the isosubstitution sites, can be useful to exclude regions for substitutions
lolimit=-10000; % Arbitrary low number
hilimit=10000;  % Arbitrary high number
dimension=2;   %x=1, y=2, z=3

% Total n of substitutions
NumTotalSubst=NumOctSubst+NumTetSubst;

ind_O1=sort([find(strcmp([atom.type],O1)) find(strcmp([atom.type],O2))]);
XYZ_labels=[atom.type]';
XYZ_data = [[atom.x]' [atom.y]' [atom.z]'];

if NumOctSubst>0;
    
    O1_atom=atom(ind_O1);
    O1_Index=1:size(O1_atom,2);%find(strcmpi(O1,strtrim(XYZ_labels(:,1))));
    O1_labels=[O1_atom.type];%XYZ_labels(O1_Index);
    O1_data=[[O1_atom.x]' [O1_atom.y]' [O1_atom.z]']; %XYZ_data(O1_Index,1:3);
    Ave_Oct_z=mean(O1_data(:,3));
    rand_O1_Index=O1_Index(randperm(length(O1_Index)));
    O1_dist_matrix = dist_matrix_atom(O1_atom,Box_dim);
    
    i=2; nOctlo=0; nOcthi=0; nOctmid=0; Oct_subst_index=[];%rand_O1_Index(1);
    while (nOctlo+nOcthi+nOctmid)<=NumOctSubst;
        ind_O2=find(strcmp([O1_atom.type],O2));
        O=intersect(ind_O2,find(O1_dist_matrix(rand_O1_Index(i),:)<minO2O2_dist));
        if length(O)<1 && nOctlo < NumOctSubst && lolimit < XYZ_data(rand_O1_Index(i),dimension) && hilimit > XYZ_data(rand_O1_Index(i),dimension);
            if nOctlo < NumOctSubst/2 && XYZ_data(ind_O1(rand_O1_Index(i)),3)<Ave_Oct_z; %&& ceil(XYZ_data(rand_O1_Index(i),3)*100)/100<=Ave_Oct_z;
                Oct_subst_index=[Oct_subst_index; rand_O1_Index(i)];
                nOctlo=nOctlo+1;
                [O1_atom(rand_O1_Index(i)).type]=O2;
            elseif nOcthi < NumOctSubst/2 && XYZ_data(ind_O1(rand_O1_Index(i)),3)>Ave_Oct_z; %&& ceil(XYZ_data(rand_O1_Index(i),3)*100)/100>=Ave_Oct_z;
                Oct_subst_index=[Oct_subst_index; rand_O1_Index(i)];
                nOcthi=nOcthi+1;
                [O1_atom(rand_O1_Index(i)).type]=O2;
            elseif (nOctlo+nOcthi+nOctmid) < NumOctSubst && XYZ_data(ind_O1(rand_O1_Index(i)),3)==Ave_Oct_z;
                Oct_subst_index=[Oct_subst_index; rand_O1_Index(i)];
                nOctmid=nOctmid+1;
                [O1_atom(rand_O1_Index(i)).type]=O2;
            end
        end
        
        if (nOctlo+nOcthi+nOctmid) == NumOctSubst
            break
        end
        if i == length(O1_data);
            break
            disp('Stopped the loop')
            pause
        end
        i=i+1;
    end
    XYZ_labels(ind_O1(Oct_subst_index))=O2;
    [atom(ind_O1(Oct_subst_index)).type]=deal(O2);
    O2_atom=atom((ind_O1(Oct_subst_index)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NumTetSubst>0;
    
    ind_T1=find(strcmp([atom.type],T1));
    ind_O2=find(strcmp([atom.type],O2));
    T1_atom=atom(ind_T1);
    assignin('caller','ind_T1',ind_T1);
    assignin('caller','T1_atom',T1_atom);
    T1O2_atom=[T1_atom O2_atom];
    T1_Index=1:size(T1_atom,2);%find(strcmpi(O1,strtrim(XYZ_labels(:,1))));
    T1_labels=[T1_atom.type];
    T1_data=[[T1_atom.x]' [T1_atom.y]' [T1_atom.z]']; %XYZ_data(O1_Index,1:3);
    Ave_Tet_z=mean(T1_data(:,3));
    rand_T1_Index=T1_Index(randperm(length(T1_Index)));
    T1_dist_matrix = dist_matrix_atom(T1_atom,Box_dim);
    T1O2_dist_matrix = dist_matrix_atom(T1_atom,O2_atom,Box_dim);% dist_matrixes_atom(T1_atom,O2_atom,Box_dim);
    
    i=2; nTetlo=0; nTethi=0; All_subst_index=Oct_subst_index; Tet_subst_index=[];%(1)= rand_T1_Index(1);
    while (nTetlo+nTethi)<=NumTetSubst;
        ind_T2=find(strcmp([T1_atom.type],T2));
        ind_T2=find(strcmp([T1_atom.type],T2));
        T=intersect(ind_T2,find(T1_dist_matrix(rand_T1_Index(i),:)<minT2T2_dist));
        TO=T1O2_dist_matrix(rand_T1_Index(i),:);
        TO=TO(TO<minT2T2_dist);
        if length(T)<1 && length(TO)<1 &&  nTetlo < NumTetSubst && lolimit < XYZ_data(rand_T1_Index(i),dimension) && hilimit > XYZ_data(rand_T1_Index(i),dimension);
            if nTetlo < NumTetSubst/2 && ceil(XYZ_data(ind_T1(rand_T1_Index(i)),3)*100)/100<=Ave_Tet_z ;
                All_subst_index=[All_subst_index; rand_T1_Index(i)];
                Tet_subst_index=[Tet_subst_index; rand_T1_Index(i)];
                nTetlo=nTetlo+1;
                [T1_atom(rand_T1_Index(i)).type]=T2;
            elseif nTethi < NumTetSubst/2 && ceil(XYZ_data(ind_T1(rand_T1_Index(i)),3)*100)/100>=Ave_Tet_z;
                All_subst_index=[All_subst_index; rand_T1_Index(i)];
                Tet_subst_index=[Tet_subst_index; rand_T1_Index(i)];
                nTethi=nTethi+1;
                [T1_atom(rand_T1_Index(i)).type]=T2;
            end
        end
        if i == length(T1_data);
            break
            disp('Stopped the loop')
            pause
        end
        i=i+1;
    end
    
    XYZ_labels(ind_T1(Tet_subst_index))=T2;
    [atom(ind_T1(Tet_subst_index)).type]=deal(T2);
    
    if nTetlo==nTethi && (nTetlo+nTethi) == NumTetSubst;
        disp('Tetrahedral substitution success!!!')
    else
        disp('Tetrahedral substitution not successful!!!')
        pause
    end
    
end

if NumOctSubst>0;
    atom_O2=atom(find(strcmpi([atom.type],O2)));
    O2_distmatrix=dist_matrix_atom(atom_O2,Box_dim);
    disp('Minimum O2O2_dist is in Å')
    min(O2_distmatrix(2:end,1))
end

if NumTetSubst>0;
    atom_T2=atom(sort([find(strcmpi([atom.type],T2)) find(strcmpi([atom.type],O2))]));
    T2_distmatrix=dist_matrix_atom(atom_T2,Box_dim);
    disp('Minimum minT2T2_dist is in Å')
    min(T2_distmatrix(2:end,1))
    
    T2O2_distmatrix=dist_matrix_atom(atom_T2,atom_O2,Box_dim);
    disp('Minimum minT2O2_dist is in Å')
    min(T2O2_distmatrix(2:end,1))
    
end

if (nOctlo==nOcthi && (nOctlo+nOcthi) == NumOctSubst) || nOctmid == NumOctSubst;
    disp('Octahedral substitution success!!!')
else
    disp('Octahedral substitution not successful!!!')
    nOctlo
    nOcthi
    nOctmid
    pause
end

%%%%%% Write structure to file %%%%%%
% write_atom_gro(atom,Box_dim,filename_out);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Plot the new structure in VMD %%%%%%
% vmd(atom(ismember([atom.type],[O1 O2 T2])),Box_dim)
%%%%%%%%%



