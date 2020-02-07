%% clayff_atom.m
% * This function tries to assign all atoms according to the clayff atom 
% types (with modified atom names by MHolmboe), with some modifications
% for edges, + ion pair potentials from Joung and Cheatham, 2008 + some 
% others...
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=clayff_atom(atom,Box_dim)
% # atom=clayff_atom(atom,Box_dim,'clayff','spc')
% # atom=clayff_atom(atom,Box_dim,'clayff','spc','edge')
%
function atom=clayff_atom(atom,Box_dim,varargin)

format compact;

if nargin >3
    ffname=varargin{1};
    watermodel=varargin{2};
else
    ffname='clayff';
    watermodel='spc/e';
end

if nargin >4
    heal=7; % with edge structure
else
    heal=1; % No healing, will run much faster
end

Heal_Al=0;
Heal_Mgo=0;
Heal_Si=0;
Heal_O=1;
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
        if strncmpi(atom(i).type,{'Si'},2);atom(i).element={'Si'};
        elseif strncmpi(atom(i).type,{'Al'},2);atom(i).element={'Al'};
        elseif strncmpi(atom(i).type,{'Mg'},2);atom(i).element={'Mg'};
        elseif strncmpi(atom(i).type,{'Fe'},2);atom(i).element={'Fe'};
        elseif strncmpi(atom(i).type,{'Ca'},2);atom(i).element={'Ca'};
        elseif strncmpi(atom(i).type,{'Ow'},2);atom(i).element={'Ow'};
        elseif strncmpi(atom(i).type,{'Hw'},2);atom(i).element={'Hw'};
        elseif strncmpi(atom(i).type,{'O'},1);atom(i).element={'O'};
        elseif strncmpi(atom(i).type,{'H'},1);atom(i).element={'H'};
        else
            [atom(i).element]=atom(i).type;
        end
    end
    
    [atom.type]=atom.element;
    [atom.fftype]=atom.element;
    
    XYZ_labels=[atom.type]';
    XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
    
    clayff_param(sort(unique([atom.type])),'SPC/E');
    
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
    
    Al=double(repmat(strncmp(XYZ_labels,'Al',1),1,length(XYZ_labels)));
    Mgo=double(repmat(strncmp(XYZ_labels,'Mgo',1),1,length(XYZ_labels)));
    Si=double(repmat(strncmp(XYZ_labels,'Si',1),1,length(XYZ_labels)));
    Cao=double(repmat(strncmp(XYZ_labels,'Cao',1),1,length(XYZ_labels)));
    O=double(repmat(strncmp(XYZ_labels,'O',1),1,length(XYZ_labels)));
    H=double(repmat(strncmp(XYZ_labels,'H',1),1,length(XYZ_labels)));
    
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
    atom = rmfield(atom,'neigh');
    for i=1:length(XYZ_labels)
        
        neigh=1;
        if strncmp(XYZ_labels(i),'Al',1)
            neigh=6;
        elseif strncmp(XYZ_labels(i),'Mgo',2)
            neigh=6;
        elseif strncmp(XYZ_labels(i),'Ca',2)
            neigh=2;
        elseif strncmp(XYZ_labels(i),'Si',1)
            neigh=4;
        elseif strncmp(XYZ_labels(i),'O',1)
            neigh=4;
        elseif strncmp(XYZ_labels(i),'H',1)
            neigh=1;
            %         elseif strncmp(XYZ_labels(i),'H',1)
            %             neigh=0;
        end
        
        k=1;j=2;
        while j < 12 && k <= neigh%atom(i).neigh=[];
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
    
    %     assignin('base','bond_matrix',bond_matrix);
    %     assignin('base','radius_limit',radius_limit);
    
    
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
                if strncmpi(atom(i).type,{'Si'},2) % Si
                    
                    Neigh_cell = sort([atom(i).neigh.type]);
                    Neighbours=strcat(Neigh_cell{:});
                    if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Si O O O O
                        atom(i).fftype={'Si'};
                    elseif length(Neigh_ind) > 4
                        disp('Si atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='Si^';
                    elseif length(Neigh_ind) == 3
                        disp('Si atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='Si_';
                        if Heal_Si == 1
                            atom(size(atom,2)+1)=atom(find(strcmp([atom.type],'O'),1));
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
                if strncmpi(atom(i).type,{'Al'},2) % Al
                    
                    if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Al O O O O O O
                        atom(i).fftype={'Al'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 5 % Al O O O O
                        atom(i).fftype={'Ale'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Al O O O O
                        atom(i).fftype={'Alt'};
                    elseif length(Neigh_ind) > 6
                        disp('Al atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='Al';
                    elseif length(Neigh_ind) > 4 && length(Neigh_ind) < 6
                        disp('Al atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='Al_';
                        if Heal_Al == 1
                            atom(size(atom,2)+1)=atom(find(strcmp([atom.type],'O'),1));
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
                if strncmpi(atom(i).type,{'Mg'},2) % Mgo
                    
                    if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Mgo O O O O O O
                        if sum(strncmpi([atom.type],'Al',2))<sum(strncmpi([atom.type],'Mg',2))
                            atom(i).fftype={'Mgh'};
                        else
                            atom(i).fftype={'Mgo'};
                        end
                    elseif length(Neigh_ind) > 6
                        disp('Mgo atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='Mgo^';
                    elseif length(Neigh_ind) > 4 && length(Neigh_ind) < 6
                        disp('Mgo atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='Mgo_';
                        if Heal_Mgo == 1
                            atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'O',1),1));
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
                if strncmpi(atom(i).type,{'Ca'},2) % Mgo
                    if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 2 % Cao O O
                        atom(i).fftype={'Cao'};
                    elseif length(Neigh_ind) > 2
                        disp('Cao atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='Cao^';
                    elseif length(Neigh_ind) < 2
                        disp('Ca atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='Ca';
                    else
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
                        i
                        %pause
                    end
                end
                
                % If the element is H
                if strncmp(atom(i).type,{'H'},1) % H
                    
                    if size(Neigh_cell,1) == 1
                        atom(i).fftype={'H'};
                    elseif length(Neigh_ind) > 1
                        disp('H atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='H^';
                        atom(i).fftype='H';
                    elseif length(Neigh_ind) < 1
                        disp('H atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='H_';
                        if Heal_H == 1
                            atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'O',1),1));
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
                if strncmp(atom(i).type,{'O'},1)
                    
                    All_Neighbours=[All_Neighbours;Neighbours];
                    if strcmp(Neighbours,'AlAlH')
                        atom(i).fftype={'Oh'};
                    elseif strcmp(Neighbours,'AlAlSi')
                        atom(i).fftype={'O'};
                    elseif strncmp(Neighbours,'HHMgo',4)
                        atom(i).fftype={'Oalhh'};
                    elseif strncmp(Neighbours,'HHMg',4)
                        atom(i).fftype={'Oalhh'};
                    elseif strcmp(Neighbours,'AlHH')
                        atom(i).fftype={'Oalhh'};
                    elseif strcmp(Neighbours,'AlHSi')
                        atom(i).fftype={'Oahs'}; % Al-OH-Si for acidic edge
                    elseif strcmp(Neighbours,'AlH')
                        atom(i).fftype={'Oalh'}; % Al-O-H or Al-O-Si
                    elseif  strncmp(Neighbours,'AlHMgo',5)
                        atom(i).fftype={'Ohmg'};
                    elseif  strncmp(Neighbours,'AlHMg',5)
                        atom(i).fftype={'Ohmg'};
                    elseif strcmp(Neighbours,'AlMgoSi')
                        atom(i).fftype={'Omg'};
                    elseif strcmp(Neighbours,'AlMgSi')
                        atom(i).fftype={'Omg'};
                    elseif strcmp(Neighbours,'HSi')
                        atom(i).fftype={'Osih'};
                    elseif strcmp(Neighbours,'AlOmg')
                        atom(i).fftype={'Odsub'};
                    elseif strcmp(Neighbours,'AlAltH')
                        atom(i).fftype={'Oalt'};
                    elseif strcmp(Neighbours,'AlAlAlt')
                        atom(i).fftype={'Oalt'};
                    elseif strcmp(Neighbours,'AlAlAl')
                        atom(i).fftype={'Oalt'};
                    elseif strcmp(Neighbours,'AltSi')
                        atom(i).fftype={'Oalt'};
                    elseif strcmp(Neighbours,'AlSi') % Al-O-H or Al-O-Si
                        if sum(ismember(find(strcmp([atom.fftype],'Alt')),[atom(i).neigh.index])) < 1
                            atom(i).fftype={'Oalsi'};
                        elseif sum(ismember(find(strcmp([atom.fftype],'Alt')),[atom(i).neigh.index])) > 0
                            atom(i).fftype={'Oalt'};
                        end
                        %                     elseif strcmp(Neighbours,'AlAl')
                        %                         atom(i).fftype={'O'};
                    elseif strcmp(Neighbours,'MgoMgoSi')
                        atom(i).fftype={'Odsub'};
                    elseif strcmp(Neighbours,'MgMgSi')
                        atom(i).fftype={'Odsub'};
                    elseif strcmp(Neighbours,'AlMgMgMg')
                        atom(i).fftype={'Oalt'};
                    elseif strcmp(Neighbours,'MgMgMgSi')
                        atom(i).fftype={'O'};
                    elseif strcmp(Neighbours,'HMgMgMg')
                        atom(i).fftype={'Oh'};
                    elseif strcmp(Neighbours,'CaCaSi')
                        atom(i).fftype={'Osi'};
                    elseif strcmp(Neighbours,'Si')
                        atom(i).fftype={'Osi'};
                    elseif strcmp(Neighbours,'SiSi')
                        atom(i).fftype={'O'};
                    elseif strcmp(Neighbours,'HH')
                        atom(i).fftype={'Ow'};
                        disp('Water?')
                    elseif length(Neigh_ind) > 2
                        disp('O atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='O^';
                    elseif length(Neigh_ind) == 1 || strcmp(Neighbours,'Si') || strcmp(Neighbours,'AlAl') 
                        disp('O atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype='O_';                       
                        if Heal_O == 1
                            try
                                atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
                            catch
                                atom(size(atom,2)+1)=atom(i);
                            end
                            atom(end).type={'H'};
                            atom(end).fftype={'H'};
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
                    if strcmp(atom(i).fftype,'Oalh') && Add_H == 1 %&& nH_extra > nH_added;
                        disp('Adding acidic H to Oalh-->Oalhh')
                        atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
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
                    if strcmp(atom(i).fftype,'Oalsi') && Add_extra_H == 1 %&& nH_extra > nH_added;
                        disp('Adding acidic H to Oalsi-->Oahs')
                        atom(i).fftype={'Oahs'};
                        atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
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

% assignin('caller','newatom',atom);

for i=1:length(unique([atom.type]))
    new_Atom_label=sort(unique([atom.type]));
    ind=ismember([atom.type],new_Atom_label(i));
    assignin('caller',strcat(char(new_Atom_label(i)),'_atom'),atom(ind));
end
ffname
watermodel

%% If the usual atom types
atom = charge_atom(atom,Box_dim,ffname,watermodel);
%% If new atom types
% atom = charge_atom(atom,Box_dim,ffname,watermodel,'more');

% Atom_label=sort(unique([atom.type]));
% clayff_param(sort(Atom_label),'SPC/E');
% no_O_label=Atom_label(~strncmp(Atom_label,'O',1));
% no_O_ind=ismember(Atom_label,no_O_label);
% atom = charge_clayff_atom(atom,Box_dim,Atom_label(no_O_ind),Charge(no_O_ind));
% atom = charge_clayff_atom(atom,Box_dim)
% Total_charge

