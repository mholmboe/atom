%% interface_atom.m
% * This function tries to assign all atoms according to the interface atom
% types (with modified atom names by MHolmboe), with some modifications for
% edges...
% * heal_iterations should be a list of numbers 1-7, corresponding to which
% assigment runs you need in order to heal Al, Mg, Si, O, H and so on...
% * The variables distance_factor and rmaxlong are related to the
% neighbour/bond cutoff radius for each atomtype
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = interface_atom(atom,Box_dim)
% * atom = interface_atom(atom,Box_dim,'interface','spc')
% * atom = interface_atom(atom,Box_dim,'interface','spc',heal_iterations)


function atom=interface_atom(atom,Box_dim,varargin)

format compact;

if nargin >3
    ffname=varargin{1};
    watermodel=varargin{2};
else
    ffname='interface';
    watermodel='tip3p';
end

if nargin>4
    heal_iterations=varargin{3};
else
    heal_iterations=1; % No healing, will run much faster
end

atom=element_atom(atom);

distance_factor=0.6;
rmaxlong=2.25;
All_Neighbours=[];

% Initialize some variables
rm_ind=[];
Heal_Al=0;
Heal_Mgo=0;
Heal_Feo=0;
Heal_Mn=0;
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
        Heal_Mgo=0
    elseif assignment_run==3
        Heal_Si=0
    elseif assignment_run==4
        Heal_O=1
    elseif assignment_run==5
        Heal_H=1
    elseif assignment_run==6
        Add_H=0
    elseif assignment_run==7
        Add_extra_H=0 % As a final step, add extra H's to a Oalhh atoms...nH_extra=16;
    end

    All_Neighbours={};
    Bond_index=zeros(4*size(size(atom,2),1),3);
    Angle_index=zeros(4*size(size(atom,2),1),4);

    [atom.element]=atom.type;

    for i=1:size(atom,2)
        if strncmpi(atom(i).type,{'Si'},2);atom(i).element={'Si'};
        elseif strncmp(atom(i).type,{'SC'},2);atom(i).element={'Si'};
        elseif strncmpi(atom(i).type,{'Alt'},3);atom(i).element={'Alt'};
        elseif strncmpi(atom(i).type,{'Al'},2);atom(i).element={'Al'};
        elseif strncmpi(atom(i).type,{'Mg'},2);atom(i).element={'Mg'};
        elseif strncmpi(atom(i).type,{'Fe'},2);atom(i).element={'Fe'};
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

    atom=bond_atom(atom,Box_dim,rmaxlong,distance_factor);
    atom=remove_sametype_bond(atom,Box_dim,Bond_index);

    assignin('caller','nBonds',nBonds);
    assignin('caller','radius_limit',radius_limit);
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','Neigh_index',Neigh_index);
    assignin('caller','dist_matrix',dist_matrix);

    i=1;nH_extra=0;
    while i <= size(atom,2)
        if mod(i,100)==1
            i-1
        end

        if strncmpi([atom(i).resname],'SOL',3)==0 && strncmpi([atom(i).resname],'ION',3)==0

            Neigh_ind=zeros(12,1);
            Neigh_vec=zeros(12,3);
            n=1;neigh=1;

            if numel([atom(i).neigh])>0

                Neigh_cell = sort([atom(i).neigh.type]);
                if length(Neigh_cell) > 0
                    Neighbours=strcat(Neigh_cell{:});
                    All_Neighbours=[All_Neighbours;{Neighbours}];
                else
                    Neighbours={'Nan'};
                end
                Neigh_ind(~any(Neigh_ind,2),:) = [];
                Neigh_vec(~any(Neigh_vec,2),:) = [];

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
                        atom(i).fftype={'Si^'};
                        atom(i).neigh.dist
                        atom(i).neigh.index
                    elseif length(Neigh_ind) == 3
                        disp('Si atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Si_'};
                        if Heal_Si == 1
                            atom(size(atom,2)+1)=atom(find(strcmp([atom.type],'O'),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+0.85*max_distance*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1.5*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1.5*Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
                    else
                        i
                        atom(i).type
                        Neighbours
                        %                 atom(i)=[];
                        %                 disp('removed atom...')
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
                            atom(size(atom,2)+1)=atom(find(strncmpi([atom.type],'O',1),1));
                            atom(end).index=size(atom,2);
                            atom(end).type={'O'};
                            atom(end).fftype={'O'};
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]-1.7*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
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
                        atom(i).fftype={'H'};
                    elseif length(Neigh_ind) > 1
                        disp('H atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'H^'};
                        atom(i).fftype={'H'};
                    elseif length(Neigh_ind) < 1
                        disp('H atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'H_'};
                        if Heal_H == 1
                            atom(size(atom,2)+1)=atom(find(strncmpi([atom.type],'O',1),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        elseif length(Neigh_ind) < 1
                            Neighbours
                            atom(i)=[];
                            disp('removed atom...')
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
                        atom(i).fftype={'Op'};
                    elseif strncmp(Neighbours,'HHMgo',4)
                        atom(i).fftype={'Oalhh'};
                    elseif strncmp(Neighbours,'HHMg',4)
                        atom(i).fftype={'Oalhh'};
                    elseif strcmp(Neighbours,'AlHH')
                        atom(i).fftype={'Oalhh'};
                    elseif strcmp(Neighbours,'AlHSi')
                        atom(i).fftype={'Oahs'}; % Al-OH-Si for acidic edge
                        %                     elseif strcmp(Neighbours,'AlAlAlH');
                        %                         atom(i).fftype={'Oahs'}; % Al-OH-Si for acidic edge
                    elseif strcmp(Neighbours,'AlH')
                        atom(i).fftype={'Oalh'}; % Al-O-H or Al-O-Si
                    elseif  strcmp(Neighbours,'AlHMgo')
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
                    elseif strcmp(Neighbours,'AlAlAl')
                        atom(i).fftype={'Op'};
                        %                         atom(i).fftype={'Oalt'};
                    elseif strcmp(Neighbours,'AlAlAlt')
                        atom(i).fftype={'Op'};
                        %                         atom(i).fftype={'Oalt'};
                    elseif strcmp(Neighbours,'AltSi')
                        atom(i).fftype={'Oalt'};
                    elseif strcmp(Neighbours,'MgSi') % New atom type for O at Si-O-Mg at the edges
                        atom(i).fftype={'Omgsi'};
                    elseif strcmp(Neighbours,'AlSi') % Al-O-H or Al-O-Si
                        if sum(ismember(find(strcmp([atom.fftype],'Alt')),[atom(i).neigh.index])) > 0
                            atom(i).fftype={'Oalt'};
                        elseif sum(ismember(find(strcmp([atom.fftype],'Alt')),[atom(i).neigh.index])) < 1
                            atom(i).fftype={'Oalsi'}; % Believe this was a special zeolite thing
                        end
                        %                 elseif strcmp(Neighbours,'AlAl')
                        %                     atom(i).fftype={'O'};
                    elseif strcmp(Neighbours,'MgoMgoSi')
                        atom(i).fftype={'Odsub'};
                    elseif strcmp(Neighbours,'MgMgSi')
                        atom(i).fftype={'Odsub'};
                    elseif strcmp(Neighbours,'MgMgMgSi') || strcmp(Neighbours,'FeoFeoFeoSi')
                        atom(i).fftype={'Op'};
                    elseif strcmp(Neighbours,'HMgMgMg') || strcmp(Neighbours,'HFeoFeoFeo')
                        atom(i).fftype={'Oh'};
                    elseif strcmp(Neighbours,'SiSi')
                        atom(i).fftype={'Ob'};
                    elseif strcmp(Neighbours,'HH')
                        atom(i).fftype={'Ow'};
                        disp('Water?')
                    elseif length(Neigh_ind) > 2
                        disp('O atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype='O^';
                        assignin('caller','atom_problem',atom);
                    elseif length(Neigh_ind) == 1 || strcmp(Neighbours,'AlAl') || strcmp(Neighbours,'AlMg')
                        if strcmp(Neighbours,'Si')
                            atom(i).fftype={'Osi'};
                        elseif strcmp(Neighbours,'Al')
                            atom(i).fftype={'Oal'};
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
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) -2*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) -2*Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
                    elseif length(Neigh_ind) == 0
                        Neighbours
                        atom(i)
                        %                     atom(i)=[];
                        disp('remove O atom...?')

                        rm_ind=[rm_ind i];

                    else
                        disp('Cannot assign this one...')
                        atom(i).type
                        i
                        Neighbours
                    end

                    %%


                    if strcmp(atom(i).fftype,'Oalh') && Add_H == 1 %&& nH_extra > nH_added;
                        disp('Adding acidic H to Oalh-->Oalhh')
                        atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
                        atom(end).index=size(atom,2);
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1/2*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1/2*Neigh.r_vec(:,3)],1)));
                        [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    end
                    %% Special thingy to add Na on the edge
                    %             if strcmp(atom(i).fftype,'Oalh') && Add_H == 1 ;%&& nH_extra > nH_added;
                    %                 disp('Adding Na!!!')
                    %                 atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
                    %                 atom(end).index=size(atom,2);
                    %                 atom(end).fftype={'Na'};
                    %                 [atom(end).type]=atom(end).fftype;
                    %                 NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]+3*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 10*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 10*Neigh.r_vec(:,3)],1)));
                    %                 [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    %             end
                    %%
                    if strcmp(atom(i).fftype,'Oalsi') && Add_extra_H == 1 ;%&& nH_extra > nH_added;
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

for i=1:length(unique([atom.type]))
    new_Atom_label=sort(unique([atom.type]));
    ind=ismember([atom.type],new_Atom_label(i));
    assignin('caller',strcat(char(new_Atom_label(i)),'_atom'),atom(ind));
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

nAtoms=size(atom,2);

ind_edgeO=find(ismember([atom.type],{'Oalhh','Oalsi','Oalhh'}));
ind_H=find(ismember([atom.type],{'H'}));
[row,col]=find(dist_matrix(ind_edgeO,:));
ind_edgeH=intersect(col,ind_H);
[atom(ind_edgeH).type]=deal({'He'});

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

for i=1:length(unique([atom.type]))
    new_Atom_label=sort(unique([atom.type]));
    ind=ismember([atom.type],new_Atom_label(i));
    assignin('caller',strcat(char(new_Atom_label(i)),'_atom'),atom(ind));
end

try
    %% If the usual atom types
    atom = charge_atom(atom,Box_dim,ffname,watermodel);
    %     pause(2)
    %     %% If new atom types
    %     atom = charge_atom(atom,Box_dim,ffname,watermodel,'adjust');
    assignin('caller','Total_charge',Total_charge);
catch
    disp('Could not set the charge...')
end

if abs(round2dec(sum([atom.charge])) - sum([atom.charge]))>0.00001
    disp('Initial total charge!')
    Total_charge=sum([atom.charge])
    [atom(1).charge]=atom(1).charge+(round2dec(sum([atom.charge])) - sum([atom.charge]));
    %     atom=tweak_charge_atom(atom,'Al');
    disp('Final total charge, modfying the first atoms charge somewhat!')
    Total_charge=sum([atom.charge])
    round2dec(sum(Total_charge),6)
    assignin('caller','Total_charge',Total_charge);
end

[C,ia,ic]=unique([atom.type],'last');
[atom(ia).type]
[atom(ia).charge]

ffname
watermodel

% assignin('caller','newatom',atom);
% assignin('caller','remove_ind',rm_ind);
assignin('caller','All_Neighbours',unique(All_Neighbours));

