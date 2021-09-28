%% substitute_atom.m
% * This scripts performs isomorphous substitution, by replacing for instance
% * some O1->O2 atomtypes and/or optionally T1->T2 atomtypes
% * varargin should be something like ,NumTetrahedralSubst,T1,T2,minT2T2_dist)
% * atom is the atom struct
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = substitute_atom(atom,Box_dim,5,'Al','Mgo',5.5) % Basic input arguments
% # atom = substitute_atom(atom,Box_dim,5,'Al','Mgo',5.5,2,'Si','Al',5.5) % Will perform both octahedral and tetrahedral replacements 
% # atom = substitute_atom(atom,Box_dim,5,'Al','Mgo',5.5,2,'Si','Al',5.5,-2.5,12.5,3) % Only subst. between z>-2.5 and z<12.5 in the z-direction (3).
%
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
% minO2O2_dist=5.5;               % Minimum O2/O2 substitutions distance in Å, you may decrease it to 4.5 if you add alot of charge
% minT2T2_dist=5.5;               % Minimum T2/T2 substitutions distance in Å, you may decrease it to 4.5 if you add alot of charge

shift_z=0;
if (sum([atom.z])/size(atom,2)) > 1
    shift_z=sum([atom.z])/size(atom,2); % Average z-position
    atom=translate_atom(atom,[0 0 -shift_z]) % Make the atom struct centrosymmetrish around z=0
end

if ~iscell(O1)
    disp('Converting O1 to cell')
    O1={O1};
end
if ~iscell(O2)
    disp('Converting O2 to cell')
    O2={O2};
end

if nargin > 6
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
    
    if ~iscell(T1)
        disp('Converting T1 to cell')
        T1={T1}
    end
    
    if ~iscell(T2)
        disp('Converting T2 to cell')
        T2={T2}
    end
else
    NumTetSubst=0;
    T1=O1; T2=O2;  % Alt   % T1 will get replaced by T2 NumTetSubst times
    minT2T2_dist=5.5;      % Minimum T2/T2 substitutions distance in Å, you may decrease it to 4.5 if you add alot of charge
end

if nargin > 10
    % Limits for the isosubstitution sites, can be useful to exclude regions for substitutions
    lolimit=varargin{5};% 35; % Arbitrary low number
    hilimit=varargin{6};% 85;  % Arbitrary high number
    dimension=varargin{7};    % meaning == 1=x, 2=y, 3=z
else
    % Limits for the isosubstitution sites, can be useful to exclude regions for substitutions
    lolimit=-1000000000;% 35; % Arbitrary low number
    hilimit=100000000;% 85;  % Arbitrary high number
    dimension=3;    % meaning == 1=x, 2=y, 3=z
end
% Total n of substitutions
NumTotalSubst=NumOctSubst+NumTetSubst;

ind_O1=sort([find(strcmp([atom.type],O1)) find(strcmp([atom.type],O2))]);
XYZ_labels=[atom.type]';
XYZ_data = [[atom.x]' [atom.y]' [atom.z]'];
O2_atom=[]; % New addition...

if NumOctSubst>0
    
    O1_atom=atom(ind_O1);
    O1_Index=1:size(O1_atom,2);%find(strcmpi(O1,strtrim(XYZ_labels(:,1))));
    O1_labels=[O1_atom.type];%XYZ_labels(O1_Index);
    O1_data=[[O1_atom.x]' [O1_atom.y]' [O1_atom.z]']; %XYZ_data(O1_Index,1:3);
    Ave_Oct_z=mean(O1_data(:,3));
    rand_O1_Index=O1_Index(randperm(length(O1_Index)));
    O1_dist_matrix = dist_matrix_atom(O1_atom,Box_dim);
    
    i=1; nOctlo=0; nOcthi=0; nOctmid=0; Oct_subst_index=[];%rand_O1_Index(1);
    while (nOctlo+nOcthi+nOctmid)<=NumOctSubst
        ind_O2=find(strcmp([O1_atom.type],O2));
        O=intersect(ind_O2,find(O1_dist_matrix(rand_O1_Index(i),:)<minO2O2_dist));
        if length(O)<1 && nOctlo<NumOctSubst && lolimit<O1_data(rand_O1_Index(i),dimension) && hilimit>O1_data(rand_O1_Index(i),dimension)
            %             O1_data(rand_O1_Index(i),2);
            if nOctlo < NumOctSubst/2 && XYZ_data(ind_O1(rand_O1_Index(i)),3)<Ave_Oct_z %&& ceil(XYZ_data(rand_O1_Index(i),3)*100)/100<=Ave_Oct_z;
                Oct_subst_index=[Oct_subst_index; rand_O1_Index(i)];
                nOctlo=nOctlo+1;
                [O1_atom(rand_O1_Index(i)).type]=O2;
            elseif nOcthi < NumOctSubst/2 && XYZ_data(ind_O1(rand_O1_Index(i)),3)>Ave_Oct_z  %&& ceil(XYZ_data(rand_O1_Index(i),3)*100)/100>=Ave_Oct_z;
                Oct_subst_index=[Oct_subst_index; rand_O1_Index(i)];
                nOcthi=nOcthi+1;
                [O1_atom(rand_O1_Index(i)).type]=O2;
            elseif (nOctlo+nOcthi+nOctmid) < NumOctSubst && XYZ_data(ind_O1(rand_O1_Index(i)),3)==Ave_Oct_z
                Oct_subst_index=[Oct_subst_index; rand_O1_Index(i)];
                nOctmid=nOctmid+1;
                [O1_atom(rand_O1_Index(i)).type]=O2;
            end
        end
        
        if (nOctlo+nOcthi+nOctmid) == NumOctSubst
            break
        end
        if i == length(O1_data)
            disp('Stopped the loop')
            %             sort(rand_O1_Index(1:i))
            nOctlo
            nOcthi
            nOctmid
            % pause(3)
            break
        end
        i=i+1;
    end
    XYZ_labels(ind_O1(Oct_subst_index))=O2;
    [atom(ind_O1(Oct_subst_index)).type]=deal(O2);
    O2_atom=atom((ind_O1(Oct_subst_index)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NumTetSubst>0

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
    if  NumOctSubst>0
        T1O2_dist_matrix = dist_matrix_atom(T1_atom,O2_atom,Box_dim);% dist_matrixes_atom(T1_atom,O2_atom,Box_dim);
    else
        T1O2_dist_matrix=[];
        Oct_subst_index=[];
    end
    i=2; nTetlo=0; nTethi=0; All_subst_index=Oct_subst_index; Tet_subst_index=[];%(1)= rand_T1_Index(1);
    while (nTetlo+nTethi)<=NumTetSubst
        ind_T2=find(strcmp([T1_atom.type],T2));
        ind_T2=find(strcmp([T1_atom.type],T2));
        T=intersect(ind_T2,find(T1_dist_matrix(rand_T1_Index(i),:)<minT2T2_dist));
        if numel(T1O2_dist_matrix)>1
            TO=T1O2_dist_matrix(rand_T1_Index(i),:);
            TO=TO(TO<minT2T2_dist);
        else
            TO=[];
        end
        if length(T)<1 && length(TO)<1 &&  nTetlo < NumTetSubst && lolimit < T1_data(rand_T1_Index(i),dimension) && hilimit > T1_data(rand_T1_Index(i),dimension)
            if nTetlo < NumTetSubst/2 && ceil(XYZ_data(ind_T1(rand_T1_Index(i)),3)*100)/100<=Ave_Tet_z
                All_subst_index=[All_subst_index; rand_T1_Index(i)];
                Tet_subst_index=[Tet_subst_index; rand_T1_Index(i)];
                nTetlo=nTetlo+1;
                [T1_atom(rand_T1_Index(i)).type]=T2;
            elseif nTethi < NumTetSubst/2 && ceil(XYZ_data(ind_T1(rand_T1_Index(i)),3)*100)/100>=Ave_Tet_z
                All_subst_index=[All_subst_index; rand_T1_Index(i)];
                Tet_subst_index=[Tet_subst_index; rand_T1_Index(i)];
                nTethi=nTethi+1;
                [T1_atom(rand_T1_Index(i)).type]=T2;
            end
        end
        if (nTetlo+nTethi) == NumTetSubst
            break
        end
        if i == length(T1_data)
            disp('Stopped the loop')
            % pause(3)
            break
        end
        i=i+1;
    end
    
    XYZ_labels(ind_T1(Tet_subst_index))=T2;
    [atom(ind_T1(Tet_subst_index)).type]=deal(T2);
    
    if nTetlo==nTethi && (nTetlo+nTethi) == NumTetSubst
        disp('Second substitution success!!!')
    else
        disp('Second substitution not optimal!!!')
        % pause(3)
    end
    
end

if abs(shift_z)>0
    atom=translate_atom(atom,[0 0 shift_z]) % Shift the atom structs z-coordinates back to the original
end

if NumOctSubst>0
    atom_O2=atom(find(strcmpi([atom.type],O2)));
    O2_distmatrix=dist_matrix_atom(atom_O2,Box_dim);
    disp('Minimum O2O2_dist is in Å')
    min(O2_distmatrix(2:end,1))
end

if NumTetSubst>0
    atom_T2=atom(sort([find(strcmpi([atom.type],T2)) find(strcmpi([atom.type],O2))]));
    T2_distmatrix=dist_matrix_atom(atom_T2,Box_dim);
    disp('Minimum minT2T2_dist is in Å')
    min(T2_distmatrix(2:end,1))
    
    if NumOctSubst>0
        T2O2_distmatrix=dist_matrix_atom(atom_T2,atom_O2,Box_dim);
        disp('Minimum minT2O2_dist is in Å')
        min(T2O2_distmatrix(2:end,1))
    end
    
end

if NumOctSubst>0
    if (nOctlo==nOcthi && (nOctlo+nOcthi) == NumOctSubst) || nOctmid == NumOctSubst
        disp('First substitution success!!!')
    else
        disp('First substitution not optimal!!!')
        try
            nOctlo
            nOcthi
            nOctmid
            assignin('caller','nOctlo',nOctlo)
            assignin('caller','nOcthi',nOcthi)
            assignin('caller','nOctmid',nOctmid)
        catch
            disp('No first subst...')
        end
        
    end
end

if NumTetSubst>0
    if (nTetlo==nTethi && (nTetlo+nTethi) == NumTetSubst)
        disp('Second substitution success!!!')
    else
        disp('Second substitution not optimal!!!')
        try
            nTetlo
            nTethi
            assignin('caller','nTetlo',nTetlo)
            assignin('caller','nTethi',nTethi)
        catch
            disp('No second subst...')
        end
        
    end
end


