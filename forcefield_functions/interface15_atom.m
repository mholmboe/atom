%% interface15_atom.m
% * This function tries to assign all atoms according to the interface15 ff
% * heal_iterations should be a list of numbers 1-7, corresponding to which
% assigment runs you need in order to heal Al, Mg, Si, O, H and so on...
% * Since the Interface FF 1.5 treats some atomtypes differently depending
% on the material (like 'CLAY_MINERALS' or 'SILICA'), you may have to pass 
% extra arguments, see the last examples below.
% * The variables distance_factor and rmaxlong are related to the
% neighbour/bond cutoff radius for each atomtype
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = interface15_atom(atom,Box_dim,'interface15','tip3p')
% # atom = interface15_atom(atom,Box_dim,'interface15','tip3p',heal_iterations)
% # atom = interface15_atom(atom,Box_dim,'interface15','tip3p',heal_iterations,'CLAY_MINERALS')
% # atom = interface15_atom(atom,Box_dim,'interface15','tip3p',heal_iterations,'SILICA')

function atom=interface15_atom(atom,Box_dim,varargin)

format compact;

if nargin>3
    ffname=varargin{1}
    watermodel=varargin{2}
else
    ffname='interface15'
    watermodel='tip3p'
end

if nargin>4
    heal_iterations=varargin{3};
else
    heal_iterations=1; % No healing, will run much faster
end

if nargin>5
    model_database=varargin{4}
else
    model_database='CLAY_MINERALS'
end

pause(2)

atom=element_atom(atom);

distance_factor=1.2;
rmaxlong=2.4;

Heal_Al=0;
Heal_Mgo=0;
Heal_Si=0;
Heal_O=0;
Heal_H=0; % I do not think this is important
Add_H=0; % Add H to get neutral edge
Add_extra_H=0; % Add H to get positive edge As a final step, add extra H's to a Oalhh atoms...nH_extra=16;
% n=7; % with edge structure
for assignment_run=heal_iterations
    
    % Heal broken bonds in the structure in steps, start from the center of the particle
    if assignment_run==1
        Heal_Al=0
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
    
    XYZ_labels=[atom.type]';
    XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
    
    All_Neighbours={};
    Bond_index=zeros(4*size(size(atom,2),1),3);
    Angle_index=zeros(4*size(size(atom,2),1),4);
    
    atom=element_atom(atom);
    
    [atom.fftype]=atom.type;
    
    XYZ_labels=[atom.type]';
    XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];
    
    Radiiproperties=load('Revised_Shannon_radii.mat');
    %     atom = rmfield(atom,'neigh');
    %     atom=bond_valence_atom(atom,Box_dim,1.25,rmax);
    %     interface15_param(sort(unique([atom.type])),'tip3p');
    
    if size(atom,2)>5000
        dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim,1.25,rmaxlong);
    else
        dist_matrix = dist_matrix_atom(atom,Box_dim,1.25,rmaxlong);
    end
    
    XYZ_radii=zeros(length(XYZ_labels),1);
    XYZ_formalcharge=zeros(length(XYZ_labels),1);
    Atom_label=sort(unique([atom.element]));
    for i=1:length(Atom_label)
        try
            ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),2));
        catch
            ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),1));
        end
        XYZ_radii(ismember([atom.element],Atom_label(i)))=median(Radiiproperties.CrysRadii(ind))';
        XYZ_formalcharge(ismember([atom.element],Atom_label(i)))=median(Radiiproperties.OxState(ind))';
    end
    
    assignin('caller','XYZ_radii',XYZ_radii);
    assignin('caller','XYZ_formalcharge',XYZ_formalcharge);
    
    XYZ_radii(XYZ_radii==0)=distance_factor;
    radius_matrix=repmat(XYZ_radii,1,length(XYZ_radii));
    radius_limit=(radius_matrix+radius_matrix')*distance_factor;
    dist_matrix(dist_matrix==0)=1000;
    bond_matrix=dist_matrix-radius_limit;
    dist_matrix(dist_matrix==1000)=0;
    bond_matrix(bond_matrix>0)=0;
    bond_matrix(bond_matrix<0)=1;
    disp('Radii+Radii limits')
    unique(XYZ_labels,'stable')
    unique(XYZ_radii,'stable')
    unique(radius_limit,'stable')
    assignin('caller','radius_limit',radius_limit);
    assignin('caller','bond_matrix',bond_matrix);
    assignin('caller','dist_matrix',dist_matrix);
    
    atom = rmfield(atom,'neigh');
    for i=1:length(XYZ_labels)
        k=1;j=1;
        bond_ind=find(bond_matrix(:,i));
        while j <= numel(bond_ind) && k <= numel(bond_ind) %<= neigh %atom(i).neigh=[];
            if bond_matrix(bond_ind(j),i)==1
                if XYZ_formalcharge(i)*XYZ_formalcharge(bond_ind(j))<0 && [atom(i).molid] == [atom(bond_ind(j)).molid]
                    [atom(i).neigh.dist(k)]=dist_matrix(bond_ind(j),i);
                    [atom(i).neigh.index(k)]=bond_ind(j);
                    [atom(i).neigh.type(k)]=XYZ_labels(bond_ind(j));
                    [atom(i).neigh.coords(k,:)]=[XYZ_data(bond_ind(j),1) XYZ_data(bond_ind(j),2) XYZ_data(bond_ind(j),3)];
                    [atom(i).neigh.r_vec(k,:)]=[X_dist(bond_ind(j),i) Y_dist(bond_ind(j),i) Z_dist(bond_ind(j),i)];
                    k=k+1;
                end
            end
            j=j+1;
        end
    end
    
    assignin('caller','Bond_index',Bond_index);
    
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
                        if strcmpi(model_database,'CLAY_MINERALS')
                            atom(i).fftype={'SY1'};
                            atom(i).type={'Si'};
                        elseif strcmpi(model_database,'SILICA')
                            Second_neighbours=0;
                            for nn=1:numel(Neigh_ind)
                                Second_neighbours=Second_neighbours+length([atom(Neigh_ind(nn)).neigh.index]);
                            end
                            if Second_neighbours==8
                                atom(i).fftype={'SC4'};
                                atom(i).type={'Si'};
                            elseif Second_neighbours==7
                                atom(i).fftype={'SC5'};
                                atom(i).type={'Si'};
                            else
                                atom(i).fftype={'SC5'};
                                atom(i).type={'Si'};
                                disp('Si has more than one undercoordinated oxygen neighbours')
                                Second_neighbours
                                i
                            end
                        end
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
                        atom(i).fftype={'AY1'};
                        atom(i).type={'Al'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 5 % Al O O O O
                        atom(i).fftype={'AYE1'};
                        atom(i).type={'Ale'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Al O O O O
                        atom(i).fftype={'AYT1'};
                        atom(i).type={'Alt'};
                    elseif length(Neigh_ind) > 6
                        disp('Al atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='Al^';
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
                        atom(i).fftype={'MY1'};
                        atom(i).type={'Mgo'};
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
                
                % If the element is H
                if strncmp(atom(i).type,{'H'},1) % H
                    
                    if size(Neigh_cell,1) == 1
                        if strcmpi(model_database,'CLAY_MINERALS')
                            atom(i).fftype={'HOY'};
                            atom(i).type={'H'};
                        elseif strcmpi(model_database,'SILICA')
                            atom(i).fftype={'HOY'};
                            atom(i).type={'H'};
                        end
                    elseif length(Neigh_ind) > 1
                        disp('H atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='H^';
                        %                         atom(i).fftype='HOY';
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
                        atom(i).fftype={'OY6'}; % Oh
                        atom(i).type={'Oh'};
                    elseif strcmp(Neighbours,'AlAlSi')
                        atom(i).fftype={'OY4'}; % Op
                        atom(i).type={'Op'};
                    elseif strncmp(Neighbours,'HHMgo',4)
                        atom(i).fftype={'OAHH'}; % Oalhh
                        atom(i).type={'Oalhh'};
                    elseif strncmp(Neighbours,'HHMg',4)
                        atom(i).fftype={'OAHH'}; % Oalhh
                        atom(i).type={'Oalhh'};
                    elseif strcmp(Neighbours,'AlHH')
                        atom(i).fftype={'OAHH'}; % Oalhh
                        atom(i).type={'Oalhh'};
                    elseif strcmp(Neighbours,'AlHSi')
                        atom(i).fftype={'OAHS'}; % Al-OH-Si for acidic edge
                        atom(i).type={'Oahsi'};
                    elseif strcmp(Neighbours,'AlH')
                        atom(i).fftype={'OAH'}; % Al-O-H or Al-O-Si
                        atom(i).type={'Oalh'};
                    elseif  strcmp(Neighbours,'AlHMgo')
                        atom(i).fftype={'OY9'}; % Ohmg
                        atom(i).type={'Ohmg'};
                    elseif  strncmp(Neighbours,'AlHMg',5)
                        atom(i).fftype={'OY9'}; % Ohmg
                        atom(i).type={'Ohmg'};
                    elseif strcmp(Neighbours,'AlMgoSi')
                        atom(i).fftype={'OY5'}; % Omg
                        atom(i).type={'Omg'};
                    elseif strcmp(Neighbours,'AlMgSi')
                        atom(i).fftype={'OY5'}; % Omg
                        atom(i).type={'Omg'};
                    elseif strcmp(Neighbours,'HSi')
                        if strcmpi(model_database,'CLAY_MINERALS')
                            atom(i).fftype={'OSH'}; % Osih
                            atom(i).type={'Osih'};
                        elseif strcmpi(model_database,'SILICA')
                            atom(i).fftype={'OC24'}; % Osih
                            atom(i).type={'Osih'};
                        end
                    elseif strcmp(Neighbours,'AlOmg')
                        atom(i).fftype={'Odsub'}; % Odsub
                        atom(i).type={'Odsub'};
                    elseif strcmp(Neighbours,'AlAltH')
                        atom(i).fftype={'OY2'}; % Oalt
                        atom(i).type={'Oalt'};
                    elseif strcmp(Neighbours,'AlAlAl')
                        atom(i).fftype={'OY2'}; % Oalt
                        atom(i).type={'Oalt'};
                    elseif strcmp(Neighbours,'AlAlAlt')
                        atom(i).fftype={'OY2'}; % Oalt
                        atom(i).type={'Oalt'};
                    elseif strcmp(Neighbours,'AltSi')
                        atom(i).fftype={'OY2'}; % Oalt
                        atom(i).type={'Oalt'};
                    elseif strcmp(Neighbours,'MgSi') % New atom type for O at Si-O-Mg at the edges
                        atom(i).fftype={'OMS'}; % Omgsi
                        atom(i).type={'Omgsi'};
                    elseif strcmp(Neighbours,'AlSi') % Al-O-H or Al-O-Si
                        if sum(ismember(find(strcmp([atom.fftype],'AYT1')),[atom(i).neigh.index])) < 1
                            atom(i).fftype={'OAS'}; % Oalsi
                            atom(i).type={'Oalsi'};
                        elseif sum(ismember(find(strcmp([atom.fftype],'AYT1')),[atom(i).neigh.index])) > 0
                            atom(i).fftype={'OY2'}; % Oalt
                            atom(i).type={'Oalt'};
                        end
%                     elseif strcmp(Neighbours,'AlAl')
%                         atom(i).fftype={'OY4'}; % Op Im guessing...
%                         atom(i).type={'Op'};
                    elseif strcmp(Neighbours,'MgoMgoSi')
                        atom(i).fftype={'Odsub'}; % Odsub
                        atom(i).type={'Odsub'};
                    elseif strcmp(Neighbours,'MgMgSi')
                        atom(i).fftype={'Odsub'}; % Odsub
                        atom(i).type={'Odsub'};
                    elseif strcmp(Neighbours,'SiSi')
                        if strcmpi(model_database,'CLAY_MINERALS')
                            atom(i).fftype={'OY1'}; % Ob
                            atom(i).type={'Ob'};
                        elseif strcmpi(model_database,'SILICA')
                            atom(i).fftype={'OC23'}; % Ob
                            atom(i).type={'Ob'};
                        end
                    elseif strcmp(Neighbours,'HH')
                        atom(i).fftype={'Ow'}; % Ow
                        disp('Water?')
                    elseif length(Neigh_ind) > 2
                        disp('O atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='O^';
                    elseif length(Neigh_ind) == 1 || strcmp(Neighbours,'AlAl') || strcmp(Neighbours,'AlMg')
                        if strcmpi(model_database,'SILICA')
                            atom(i).fftype={'OC25'}; % Osi
                            atom(i).type={'Osi'};
                            atom(Neigh_ind).fftype={'SC5'};
                            atom(Neigh_ind).type={'Si'};
                        else
                            disp('O atom under coordinated')
                            i
                            Neighbours
                            atom(i).fftype='O_';
                        end
                        if Heal_O == 1
                            try
                                 atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
                            catch
                                atom(size(atom,2)+1)=atom(i);
                            end
                            atom(end).type={'H'};
                            atom(end).fftype={'H'};
                            atom(end).element={'H'};
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
                        %                 pause
                        
                    else
                        disp('Cannot assign this one...')
                        atom(i).type
                        i
                        Neighbours
                    end
                    
                    if strcmp(atom(i).fftype,'OAH') && Add_H == 1 %&& nH_extra > nH_added;
                        disp('Adding acidic H to Oalh-->Oalhh')
                        atom(size(atom,2)+1)=atom(find(strncmp([atom.fftype],'HOY',1),1));
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
                    %                 NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]+3*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 10*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 10*Neigh.r_vec(:,3)],1)));
                    %                 [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    %             end
                    %
                    if strcmp(atom(i).fftype,'OAS') && Add_extra_H == 1 %&& nH_extra > nH_added;
                        disp('Adding acidic H to Oalsi-->Oahs')
                        atom(i).fftype={'OAHS'};
                        atom(size(atom,2)+1)=atom(find(strncmp([atom.fftype],'HOY',1),1));
                        atom(end).index=size(atom,2);
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 0*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 0*Neigh.r_vec(:,3)],1)));
                        [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    end
                end
            end
        end
        i=i+1;
        %[atom.type]=atom.fftype;
    end
end

[Y,I]=sort(Bond_index(:,1));
Bond_index=Bond_index(I,:);
Bond_index = unique(Bond_index,'rows','stable');

[Y,I]=sort(Angle_index(:,2));
Angle_index=Angle_index(I,:);
Angle_index = unique(Angle_index,'rows','stable');

Bond_index(~any(Bond_index,2),:) = [];
Angle_index(~any(Angle_index,2),:) = [];

for i=1:length(unique([atom.fftype]))
    new_Atom_label=sort(unique([atom.fftype]));
    ind=ismember([atom.fftype],new_Atom_label(i));
    assignin('caller',strcat(char(new_Atom_label(i)),'_atom'),atom(ind));
end

nAtoms=size(atom,2);

ind_edgeO=find(ismember([atom.fftype],{'OAHH','OAS','OAHH'}));
ind_H=find(ismember([atom.fftype],{'HOY'}));
[row,col]=find(bond_matrix(ind_edgeO,:));
ind_edgeH=intersect(col,ind_H);
[atom(ind_edgeH).fftype]=deal({'HE'});

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

for i=1:length(unique([atom.type]))
    new_Atom_label=sort(unique([atom.type]));
    ind=ismember([atom.type],new_Atom_label(i));
    assignin('caller',strcat(char(new_Atom_label(i)),'_atom'),atom(ind));
end

% try
%     %% If the usual atom types
%     atom = charge_atom(atom,Box_dim,ffname,watermodel);
%     %% If new atom types
%     atom_alternative = charge_atom(atom,Box_dim,ffname,watermodel,'adjust');
%     assignin('caller','atom_alternative',atom_alternative);
% catch
%     disp('Could not set the charge...')
% end
atom = check_interface15_charge(atom,model_database);

