%% clayff_2004_atom.m
% * This function tries to assign all atoms according to the clayff atom types (with modified atom names by MHolmboe), with some modifications for edges...
% * This clayffmod function was modified accoring to Clays and Clay Minerals, Vol. 64, No. 4, 452?471, 2016. STRUCTURE
%
%% Version
% 2.03
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=clayff_atom(atom,Box_dim)
% # atom=clayff_atom(atom,Box_dim,'clayff','spc')
% # atom=clayff_atom(atom,Box_dim,'clayff','spc','edge')
%
function atom=clayff_2004_atom(atom,Box_dim,varargin)
 
format compact;

atom=element_atom(atom);

if nargin >3
    ffname=varargin{1};
    watermodel=varargin{2};
else
    ffname='clayff_2004';
    watermodel='spc/e';
end

if nargin >4
    heal=7; % with edge structure
else
    heal=1; % No healing, will run much faster
end

% Initialize some variables
Heal_Al=0;
Heal_Mgo=0;
Heal_Si=0;
Heal_O=0;
Heal_H=0; % I do not think this is important
Add_H=0; % Add H to get neutral edge
Add_extra_H=0; % Add H to get positive edge As a final step, add extra H's to a Oalhh atoms...nH_extra=16;
% n=7; % with edge structure
for assignment_run=1:heal
    
    % Heal broken bonds in the structure in steps, start from the center of the particle
    if assignment_run==1
        Heal_Al=1
    elseif assignment_run==2
        Heal_Mgo=1
    elseif assignment_run==3
        Heal_Si=1
    elseif assignment_run==4
        Heal_O=1
    elseif assignment_run==5
        Heal_H=1
    elseif assignment_run==6
        Add_H=1
    elseif assignment_run==7
        Add_extra_H=0 % As a final step, add extra H's to a Oalhh atoms...nH_extra=16;
    end
    
    All_Neighbours={};
    Bond_index=zeros(4*size(size(atom,2),1),3);
    Angle_index=zeros(4*size(size(atom,2),1),4);
    
    [atom.element]=atom.type;
    
    for i=1:size(atom,2)
        if strncmpi(atom(i).type,{'Si'},2);atom(i).element={'st'};
        elseif strncmpi(atom(i).type,{'Al'},2);atom(i).element={'ao'};
        elseif strncmpi(atom(i).type,{'Mg'},2);atom(i).element={'mgo'};
        elseif strncmpi(atom(i).type,{'Fe'},2);atom(i).element={'feo'};
        elseif strncmpi(atom(i).type,{'Ow'},2);atom(i).element={'Ow'};
        elseif strncmpi(atom(i).type,{'Hw'},2);atom(i).element={'Hw'};
        elseif strncmpi(atom(i).type,{'O'},1);atom(i).element={'o'};
        elseif strncmpi(atom(i).type,{'H'},1);atom(i).element={'ho'};
        else
            [atom(i).element]=atom(i).type;
        end
    end
    
    [atom.type]=atom.element;
    [atom.fftype]=atom.element;
    
    XYZ_labels=[atom.type]';
    XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
    
    clayff_2004_param(sort(unique([atom.type])),'SPC/E');
    
    dist_matrix = dist_matrix_atom(atom,Box_dim);
    
    XYZ_num=zeros(length(XYZ_labels),1);
    Atom_label=sort(unique([atom.type]));
    for i=1:length(Atom_label)
        XYZ_num(ismember([atom.type],Atom_label(i)))=i;
    end
    
    XYZ_sigma=Sigma(XYZ_num(:))';
    XYZ_epsilon=Epsilon(XYZ_num(:))';
    XYZ_mass=Masses(XYZ_num(:))';
    XYZ_charge=Charge(XYZ_num(:))';
    
    distance=0.4;
    XYZ_sigma(XYZ_sigma==0)=distance;
    radius_matrix=repmat(XYZ_sigma,1,length(XYZ_sigma));
    radius_limit=(radius_matrix+radius_matrix')*distance;
    
    Al=double(repmat(strncmp(XYZ_labels,'at',1),1,length(XYZ_labels)));
    Mgo=double(repmat(strncmp(XYZ_labels,'mgo',1),1,length(XYZ_labels)));
    Si=double(repmat(strncmp(XYZ_labels,'st',1),1,length(XYZ_labels)));
    Cao=double(repmat(strncmp(XYZ_labels,'cao',1),1,length(XYZ_labels)));
    O=double(repmat(strncmp(XYZ_labels,'o',1),1,length(XYZ_labels)));
    H=double(repmat(strncmp(XYZ_labels,'ho',1),1,length(XYZ_labels)));
    
    AlO=(Al+O'>1)+(Al'+O>1);     radius_limit(AlO>0)=2.8;
    CaoO=(Cao+O'>1)+(Cao'+O>1);  radius_limit(CaoO>0)=2.8;
    MgoO=(Mgo+O'>1)+(Mgo'+O>1);  radius_limit(MgoO>0)=2.8;
    SiO=(Si+O'>1)+(Si'+O>1);     radius_limit(SiO>0)=2.1;
    HO=(H+O'>1)+(H'+O>1);        radius_limit(HO>0)=1.25;
    
    % To changes the radius manually, do it like this!!
    % radius_limit(261,256)=2.3;
    % radius_limit(251,256)=2.3;
    % radius_limit(272,278)=3;
    
    bond_matrix=(dist_matrix-radius_limit)<0;
    bond_matrix(logical(eye(size(bond_matrix)))) = 0;
    
    [bond_r,bond_ind]=sort(dist_matrix(:,:));
    bond_labels=XYZ_labels(bond_ind);
    try
    atom = rmfield(atom,'neigh');
    catch
        disp('no atom.neigh')
    end
    for i=1:length(XYZ_labels)
        
        if strncmp(XYZ_labels(i),'ao',1)
            neigh=6;
        elseif strncmp(XYZ_labels(i),'mgo',2)
            neigh=6;
        elseif strncmp(XYZ_labels(i),'cao',3)
            neigh=2;
        elseif strncmp(XYZ_labels(i),'st',1)
            neigh=4;
        elseif strncmp(XYZ_labels(i),'o',1)
            neigh=4;
        elseif strncmp(XYZ_labels(i),'ho',1)
            neigh=1;
%         elseif strncmp(XYZ_labels(i),'h',1)
%             neigh=0;
        end
        
        k=1;j=2;
        while j < 12 && k <= neigh %atom(i).neigh=[];
            if bond_matrix(bond_ind(j,i),i)==1
                if XYZ_charge(i)*XYZ_charge(bond_ind(j,i))<0
                    [atom(i).neigh.dist(k)]=bond_r(j,i);
                    [atom(i).neigh.index(k)]=bond_ind(j,i);
                    [atom(i).neigh.type(k)]=bond_labels(j,i);
                    [atom(i).neigh.coords(k,:)]=[XYZ_data(bond_ind(j,i),1) XYZ_data(bond_ind(j,i),2) XYZ_data(bond_ind(j,i),3)];
                    [atom(i).neigh.r_vec(k,:)]=[X_dist(bond_ind(1,i),bond_ind(j,i)) Y_dist(bond_ind(1,i),bond_ind(j,i)) Z_dist(bond_ind(1,i),bond_ind(j,i))];
                    %[atom(i).neigh.r_vec(k,:)]=[XYZ_data(bond_ind(1,i),1) XYZ_data(bond_ind(1,i),2) XYZ_data(bond_ind(1,i),3)]...
                    %    - [XYZ_data(bond_ind(j,i),1) XYZ_data(bond_ind(j,i),2) XYZ_data(bond_ind(j,i),3)];
                    k=k+1;
                end
            end
            j=j+1;
        end
    end
    
    b=1;a=1;i=1;nH_extra=0;
    while i <= size(atom,2)
        if mod(i,100)==1
            i-1
        end
        
        if strncmpi([atom(i).resname],'SOL',3)==0 && strncmpi([atom(i).resname],'ION',3)==0
            Neigh_ind=zeros(12,1);
            Neigh_vec=zeros(12,3);
            n=1;neigh=1;
            
            if numel([atom(i).neigh])>0
                index=[atom(i).neigh.index];
                if sum(index)>0
                    for k=1:length(index)
                        j=atom(i).neigh.index(k);
                        r=atom(i).neigh.dist(k);
                        
                        if i < size(radius_limit,2)
                            max_distance=radius_limit(j,i);
                        else
                            disp('Manually setting max_distance to...')
                            max_distance=1.8
                        end
                        if r > 0 && r < max_distance/3
                            r
                            disp('Atoms too close')
                            format compact;
                            format short;
                            i
                            j
                            atom(i).type
                            atom(j).type
                            [atom(i).x atom(i).y atom(i).z]
                            [atom(j).x atom(j).y atom(j).z]
                        end
                        
                        if r > max_distance/3 && r < max_distance
                            
                            neigh=neigh+1;
                            
                            if i < j
                                Bond_index(b,1)= i;
                                Bond_index(b,2)= j;
                            else
                                Bond_index(b,1)= j;
                                Bond_index(b,2)= i;
                            end
                            
                            Bond_index(b,3)= r;
                            Neigh_ind(n,1)= j;
                            Neigh_vec(n,1:3)= [atom(i).neigh.r_vec(k,1) atom(i).neigh.r_vec(k,2) atom(i).neigh.r_vec(k,3)];
                            % Neigh_vec(n,1:3)= [atom(i).neigh(k).r_vec(1) atom(i).neigh(k).r_vec(2) atom(i).neigh(k).r_vec(3)];
                            b=b+1;
                            n=n+1;
                        end
                    end
                end
                
                Neigh_cell = sort([atom(i).neigh.type]);
                if length(Neigh_cell) > 0; Neighbours=strcat(Neigh_cell{:});
                else Neighbours={'Nan'}; end
                Neigh_ind(~any(Neigh_ind,2),:) = [];
                Neigh_vec(~any(Neigh_vec,2),:) = [];
                
                size(Neigh_ind,1);
                for v=1:size(Neigh_ind,1)
                    for w=1:size(Neigh_ind,1)
                        angle=rad2deg(atan2(norm(cross(Neigh_vec(v,:),Neigh_vec(w,:))),dot(Neigh_vec(v,:),Neigh_vec(w,:))));
                        if angle > 0 && angle < 160
                            if v < w
                                Angle_index(a,1)= Neigh_ind(v,1);
                                Angle_index(a,2)= i;
                                Angle_index(a,3)= Neigh_ind(w,1);
                                Angle_index(a,4)= angle;
                                a=a+1;
                            else
                                Angle_index(a,1)= Neigh_ind(w,1);
                                Angle_index(a,2)= i;
                                Angle_index(a,3)= Neigh_ind(v,1);
                                Angle_index(a,4)= angle;
                                a=a+1;
                            end
                        end
                    end
                end
                
                Bond_index=sortrows(Bond_index);
                Angle_index=sortrows(Angle_index);
                
                % If the element is Si
                if strncmpi(atom(i).type,{'st'},2) % Si
                    
                    Neigh_cell = sort([atom(i).neigh.type]);
                    Neighbours=strcat(Neigh_cell{:});
                    if sum(strncmp({'o'},[atom(i).neigh.type],1)) == 4 % Si O O O O
                        atom(i).fftype={'st'};
                    elseif length(Neigh_ind) > 4
                        disp('Si atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='st^';
                    elseif length(Neigh_ind) == 3
                        disp('Si atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='st_';
                        if Heal_Si == 1
                            atom(size(atom,2)+1)=atom(find(strcmp([atom.type],'o'),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+0.75*max_distance*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1.5*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1.5*Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
                    else
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
                        i
                    end
                    %             else
                    %                 Neighbours
                    %                 atom(i)=[];
                    %                 disp('removed atom...')
                    %                 i
                    %             end
                end
                
                % If the element is Al
                if strncmpi(atom(i).type,{'ao'},2) % Al
                    
                    if sum(strncmp({'o'},[atom(i).neigh.type],1)) == 6 % Al O O O O O O
                        atom(i).fftype={'ao'};
                    elseif sum(strncmp({'o'},[atom(i).neigh.type],1)) == 5 % Al O O O O
                        atom(i).fftype={'ae'};
                    elseif sum(strncmp({'o'},[atom(i).neigh.type],1)) == 4 % Al O O O O
                        atom(i).fftype={'at'};
                    elseif length(Neigh_ind) > 6
                        disp('Al atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='ao';
                    elseif length(Neigh_ind) > 4 && length(Neigh_ind) < 6
                        disp('Al atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='ae_';
                        if Heal_Al == 1
                            atom(size(atom,2)+1)=atom(find(strcmp([atom.type],'o'),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+0.75*max_distance*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
                    else
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
                        i
                        
                        
                    end
                    %             else
                    %                 Neighbours
                    %                 atom(i)=[];
                    %                 disp('removed atom...')
                    %                 i
                    %             end
                    
                end
                
                % If the element is Mg
                if strncmpi(atom(i).type,{'mg'},2) % Mgo
                    
                    if sum(strncmp({'o'},[atom(i).neigh.type],1)) == 6 % Mgo O O O O O O
                        if sum(strncmpi([atom.type],'ao',2))<sum(strncmpi([atom.type],'mg',2))
                            atom(i).fftype={'mgh'};
                        else
                            atom(i).fftype={'mgo'};
                        end
                    elseif length(Neigh_ind) > 6
                        disp('mgo atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='mgo^';
                    elseif length(Neigh_ind) > 4 && length(Neigh_ind) < 6
                        disp('mgo atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='mgo_';
                        if Heal_Mgo == 1
                            atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'o',1),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+0.75*max_distance*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
                    else
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
                        i
                        %pause
                    end
                    
                end
                
                % If the element is Ca
                if strncmpi(atom(i).type,{'cao'},3) % cao
                    if sum(strncmp({'o'},[atom(i).neigh.type],1)) == 2 % Cao O O
                        atom(i).fftype={'cao'};
                    elseif length(Neigh_ind) > 2
                        disp('cao atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='Cao^';
                    elseif length(Neigh_ind) < 2
                        disp('ca atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='ca';
                    else
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
                        i
                        %pause
                    end
                end
                
                % If the element is H
                if strncmp(atom(i).type,{'ho'},1) % H
                    
                    if size(Neigh_cell,1) == 1
                        atom(i).fftype={'ho'};
                    elseif length(Neigh_ind) > 1
                        disp('H atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='h^';
                        atom(i).fftype='h';
                    elseif length(Neigh_ind) < 1
                        disp('H atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='h_';
                        if Heal_H == 1
                            atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'o',1),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        elseif length(Neigh_ind) < 1
                            Neighbours
                            %                     atom(i)=[];
                            %                     disp('removed atom...')
                            i
                            %pause
                        end
                        
                    end
                    
                end
                
                % If the element is O
                if strncmp(atom(i).type,{'o'},1)
                    
                    All_Neighbours=[All_Neighbours;Neighbours];                
                    if strncmp(Neighbours,'aoaoho',5)
                        atom(i).fftype={'oh'};
                    elseif strncmp(Neighbours,'aoaost',5)
                        atom(i).fftype={'op'};
                    elseif strncmp(Neighbours,'hohomgo',6)
                        atom(i).fftype={'oahh'};
                    elseif strncmp(Neighbours,'hohomg',6)
                        atom(i).fftype={'oahh'};
                    elseif strncmp(Neighbours,'aohoho',6)
                        atom(i).fftype={'oahh'};
                    elseif strncmp(Neighbours,'aohost',6)
                        atom(i).fftype={'oahs'}; % Al-OH-Si for acidic edge
                    elseif  strncmp(Neighbours,'aohomg',5) % 'aohmgo'
                        atom(i).fftype={'ohs'};
                    elseif strncmp(Neighbours,'aomgst',5)
                        atom(i).fftype={'obts'};
                    elseif strncmp(Neighbours,'aoho',4)
                        atom(i).fftype={'oah'}; % Al-O-H or Al-O-Si
                    elseif strncmp(Neighbours,'aomgost',6)
                        atom(i).fftype={'obos'};
                    elseif strncmp(Neighbours,'host',4)
                        atom(i).fftype={'ohst'};
                    elseif strncmp(Neighbours,'aoomg',5)
                        atom(i).fftype={'obss'};
                    elseif strncmp(Neighbours,'aoatho',6)
                        atom(i).fftype={'obts'};
                    elseif strncmp(Neighbours,'aoaoat',6)
                        atom(i).fftype={'obts'};
                    elseif strncmp(Neighbours,'aoat',4)
                        atom(i).fftype={'obts'};
                    elseif strncmp(Neighbours,'atst',4)
                        atom(i).fftype={'oat'};
                    elseif strncmp(Neighbours,'aost',4) % Al-O-H or Al-O-Si
                        if sum(ismember(find(strcmp([atom.fftype],'at')),[atom(i).neigh.index])) < 1
                            atom(i).fftype={'oas'};
                        elseif sum(ismember(find(strcmp([atom.fftype],'at')),[atom(i).neigh.index])) > 0
                            atom(i).fftype={'oat'};
                        end
                    elseif strncmp(Neighbours,'aoao',4)
                        atom(i).fftype={'o'};
                    elseif strncmp(Neighbours,'mgomgost',7)
                        atom(i).fftype={'odsub'};
                    elseif strncmp(Neighbours,'stst',4)
                        atom(i).fftype={'ob'}; % Old version atom(i).fftype={'O'};
                    elseif strncmp(Neighbours,'hoho',4)
                        atom(i).fftype={'o*'};
                        disp('Water?')
                    elseif length(Neigh_ind) > 2
                        disp('O atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='o^';
                    elseif length(Neigh_ind) == 1 || strncmp(Neighbours,'aoao',4) 
                        disp('O atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='o_';
                        if Heal_O == 1
                            
                            try
                                atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'ho',1),1));
                            catch
                                atom(size(atom,2)+1)=atom(i);
                            end
                            atom(end).type={'ho'};
                            atom(end).fftype={'ho'};
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) -2*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) -2*Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
                    elseif length(Neigh_ind) == 0
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
                        i
                        %pause
                    end
                    if strcmp(atom(i).fftype,'oalh') && Add_H == 1 %&& nH_extra > nH_added;
                        disp('Adding acidic H to Oalh-->Oalhh')
                        atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'h',1),1));
                        atom(end).index=size(atom,2);
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1/2*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1/2*Neigh.r_vec(:,3)],1)));
                        [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    end
                    % Special thingy to add Na on the edge
                    %             if strcmp(atom(i).fftype,'Oalh') && Add_H == 1 ;%&& nH_extra > nH_added;
                    %                 disp('Adding Na!!!')
                    %                 atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
                    %                 atom(end).index=size(atom,2);
                    %                 atom(end).fftype={'Na'};
                    %                 [atom(end).type]=atom(end).fftype;
                    %                 Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                    %                 NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]+3*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 10*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 10*Neigh.r_vec(:,3)],1)));
                    %                 [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    %             end
                    %%
                    if strcmp(atom(i).fftype,'oalsi') && Add_extra_H == 1 %&& nH_extra > nH_added;
                        disp('Adding acidic h to oalsi-->oahs')
                        atom(i).fftype={'oahs'};
                        atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'h',1),1));
                        atom(end).index=size(atom,2);
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 0*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 0*Neigh.r_vec(:,3)],1)));
                        [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    end
                end
            end
        end
        i=i+1;
        [atom.type]=atom.fftype;
    end
end

if sum(strncmp([atom.type],'Ow',2))>0
    ind_Ow=find(strncmpi([atom.type],'Ow',2));
    ind_Hw=sort([ind_Ow+1 ind_Ow+2]);
    [atom(ind_Hw).type]=deal({'Hw'});
elseif sum(strncmp([atom.type],'OW',2))>0
    ind_Ow=find(strncmpi([atom.type],'OW',2));
    ind_Hw1=sort(ind_Ow+1);
    ind_Hw2=sort(ind_Ow+2);
    [atom(ind_Hw1).type]=deal({'HW1'});
    [atom(ind_Hw2).type]=deal({'HW2'});
end

[Y,I]=sort(Bond_index(:,1));
Bond_index=Bond_index(I,:);
Bond_index = unique(Bond_index,'rows','stable');

[Y,I]=sort(Angle_index(:,2));
Angle_index=Angle_index(I,:);
Angle_index = unique(Angle_index,'rows','stable');

Bond_index(~any(Bond_index,2),:) = [];
Angle_index(~any(Angle_index,2),:) = [];

assignin('caller','newatom',atom);

for i=1:length(unique([atom.type]))
    new_Atom_label=sort(unique([atom.type]));
    ind=ismember([atom.type],new_Atom_label(i));
    assignin('caller',strcat(char(new_Atom_label(i)),'_atom'),atom(ind));
end
ffname
watermodel
atom = charge_atom(atom,Box_dim,ffname,watermodel);

% Atom_label=sort(unique([atom.type]));
% clayff_param(sort(Atom_label),'SPC/E');
% no_O_label=Atom_label(~strncmp(Atom_label,'O',1));
% no_O_ind=ismember(Atom_label,no_O_label);
% atom = charge_clayff_atom(atom,Box_dim,Atom_label(no_O_ind),Charge(no_O_ind));
% atom = charge_clayff_atom(atom,Box_dim)
% Total_charge

