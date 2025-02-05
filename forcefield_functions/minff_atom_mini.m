%% minff_atom.m
% * This function tries to assign all atoms according to the MINFF atom
% types. The variables distance_factor and rmaxlong are related to the
% neighbour/bond cutoff radius for each atomtype
%
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=minff_atom(atom,Box_dim)

function atom=minff_atom_mini(atom,Box_dim,varargin)
%%
format compact;

distance_factor=0.6;
rmaxlong=2.45;

atom=element_atom(atom);

% Initialize some variables
rm_ind=[];
for assignment_run=1:2

    All_Neighbours=[];
    [atom.element]=atom.type;

    for i=1:size(atom,2)
        if strncmpi(atom(i).type,{'Si'},2);atom(i).element={'Si'};
        elseif strncmp(atom(i).type,{'SC'},2);atom(i).element={'Si'};
        elseif strncmpi(atom(i).type,{'Ale'},3);atom(i).element={'Ale'};
        elseif strncmpi(atom(i).type,{'Alt'},3);atom(i).element={'Alt'};
        elseif strncmpi(atom(i).type,{'Al'},2);atom(i).element={'Al'};
        elseif strncmpi(atom(i).type,{'Mg'},2);atom(i).element={'Mg'};
        elseif strncmpi(atom(i).type,{'Fee'},3);atom(i).element={'Fee'};
        elseif strncmpi(atom(i).type,{'Fet'},3);atom(i).element={'Fet'};
        elseif strncmpi(atom(i).type,{'Fe'},2);atom(i).element={'Fe'};
        elseif strncmpi(atom(i).type,{'F'},1);atom(i).element={'F'};
        elseif strncmpi(atom(i).type,{'Li'},2);atom(i).element={'Li'};
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

    atom=bond_atom(atom,Box_dim,rmaxlong,distance_factor); % Genereates the variable below
    assignin('caller','nBonds',nBonds);
    assignin('caller','radius_limit',radius_limit);
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','Neigh_index',Neigh_index);
    assignin('caller','dist_matrix',dist_matrix);

    i=1;
    while i <= size(atom,2)

            nNeigh=numel([atom(i).neigh.index]);

            if nNeigh>0

                Neigh_cell = sort([atom(i).neigh.type]);
                if length(Neigh_cell) > 0
                    Neighbours=strcat(Neigh_cell{:});
                else
                    Neighbours={'Nan'};
                end

                % If the element is Li
                if strncmpi(atom(i).type,{'Li'},2) % Lio

                    if numel([atom(i).neigh.type]) == 6 %% sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Lio O O O O O O

                        atom(i).fftype={'Lio'};

                    elseif numel([atom(i).neigh.type]) == 4 %% sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Lio O O O O O O

                        atom(i).fftype={'Lio'};

                    elseif nNeigh > 6
                        disp('Lio atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Lio_ov'};
                        rm_ind=[rm_ind max(atom(i).neigh.index)];
                        pause(1) %
                    elseif nNeigh > 4 && nNeigh < 6
                        disp('Lio atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Lio_un'};

                    else
                        Neighbours
                        i
                    end

                end


                % If the element is Si
                if strncmpi(atom(i).type,{'Si'},2) % Si
                    if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Si O O O O
                        atom(i).fftype={'Sit'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 3 % Si O O O
                        atom(i).fftype={'Site'}; % As in Stishovite
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Si O O O O O O
                        atom(i).fftype={'Sio'}; % As in Stishovite
                    elseif nNeigh > 4
                        disp('Si atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Si_ov'};
                        atom(i).neigh.dist
                        atom(i).neigh.index

                        rm_ind=[rm_ind max(atom(i).neigh.index)];
                        pause(1) %
                    elseif nNeigh < 4 % == 3
                        disp('Si atom under coordinated')
                        i
                        nNeigh
                        Neighbours
                        atom(i).fftype={'Si_un'};
                    else
                        i
                        Neighbours
                    end
                end

                % If the element is Al
                if strncmpi(atom(i).type,{'Al'},2) % Al

                    if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Al O O O O O O
                        atom(i).fftype={'Al'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 5 % Al O O O O
                        i
                        atom(i).fftype={'Ale'}
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Al O O O O
                        atom(i).fftype={'Alt'};
                    elseif nNeigh > 6
                        disp('Al atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Al'};
                        rm_ind=[rm_ind max(atom(i).neigh.index)];
                        pause(1) %
                    elseif nNeigh >= 4 && nNeigh < 6
                        disp('Al atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Al_un'};
                    else
                        i
                        Neighbours
                    end
                end

                % If the element is Mg
                if strncmpi(atom(i).type,{'Mg'},2) % Mgo

                    if numel([atom(i).neigh.type]) == 6 %% sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Mgo O O O O O O
                        if sum(strncmpi([atom.type],'Mg',2))>sum(strncmpi([atom.type],'Si',2))
                            atom(i).fftype={'Mgo'}; % Ex Fosterite
                            if sum(strncmpi([atom.type],'Mg',2))<sum(strncmpi([atom.type],'H',1))
                                atom(i).fftype={'Mgh'}; % Ex. Brucite
                            end
                        elseif sum(strncmpi([atom.type],'Mg',2))<sum(strncmpi([atom.type],'Si',2))
                            atom(i).fftype={'Mgh'}; % Talc, Hectorite
                            if sum(strncmpi([atom.type],'Al',2))>sum(strncmpi([atom.type],'Mg',2))
                                atom(i).fftype={'Mgo'}; % Mica, smectite
                            end
                        end
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 5
                        atom(i).fftype={'Mge'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4
                        atom(i).fftype={'Mgo'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 2
                        atom(i).fftype={'Mgh'};
                    elseif nNeigh > 6
                        disp('Mgo atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Mgo_ov'};
                        rm_ind=[rm_ind max(atom(i).neigh.index)];
                        pause(1) %
                    elseif nNeigh > 4 && nNeigh < 6
                        disp('Mgo atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Mgo_un'};
                        atom(i).fftype={'Mgo'};
                    else
                        Neighbours
                        i
                    end

                end


                % If the element is Ca
                if strncmpi(atom(i).type,{'Ca'},2) % Cao

                    if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Cao O O O O O O
                        atom(i).fftype={'Cao'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4
                        atom(i).fftype={'Cah'};
                        %                     elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 2
                        %                         atom(i).fftype={'Cah'};
                    elseif sum(strncmp({'F'},[atom(i).neigh.type],1)) == 8
                        disp('Ca in CaF2 Fluorite?')
                        atom(i).fftype={'Cah'};
                        %                     elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 2
                        %                         atom(i).fftype={'Cah'};
                    elseif nNeigh > 6
                        disp('Ca atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Cao_ov'};
                        rm_ind=[rm_ind max(atom(i).neigh.index)];
                        pause(1) %
                    else
                        Neighbours
                    end

                end


                % If the element is Fe
                if strncmpi(atom(i).type,{'Fe'},2) % Fe
                    if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Fe O O O O O O
                        if sum([atom(i).neigh.dist])<12.7
                            atom(i).fftype={'Feo'};
                        else
                            atom(i).fftype={'Fe2'};
                        end
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 5 % Fe O O O O
                        atom(i).fftype={'Fee'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Fe O O O O
                        if sum([atom(i).neigh.dist])<7.5
                            atom(i).fftype={'Fet'};
                        else
                            atom(i).fftype={'Fe2'};
                        end
                    elseif nNeigh > 6
                        disp('Fe atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Fe_ov'};
                        rm_ind=[rm_ind max(atom(i).neigh.index)];
                        pause(1) %
                    elseif nNeigh > 4 && nNeigh < 6
                        disp('Fe atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Fe_un'};
                    else
                        Neighbours
                        i
                    end

                elseif strcmpi(atom(i).type,{'F'}) > 0% Fs

                    if numel([atom(i).neigh.type]) == 3

                        atom(i).fftype={'Fs'};

                    elseif nNeigh == 4
                        disp('Fs atom as in CaF2 - Fluorite?')
                        atom(i).fftype={'Fs'};

                    elseif nNeigh > 4
                        disp('Fso atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Fs_ov'};
                        rm_ind=[rm_ind max(atom(i).neigh.index)];
                        pause(1) %
                    elseif nNeigh > 0 && nNeigh < 3
                        disp('Fso atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Fs_un'};

                    else
                        Neighbours
                        i
                    end

                end

                % If the element is Ti
                if strncmpi(atom(i).type,{'Ti'},2) % Ti
                    %                     Neigh_cell = sort([atom(i).neigh.type]);
                    %                     Neighbours=strcat(Neigh_cell{:});
                    if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Ti O O O O
                        atom(i).fftype={'Tit'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Ti O O O O
                        atom(i).fftype={'Tio'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) > 6 % Ti O O O O
                        disp('Ti atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Ti_ov'};
                        atom(i).neigh.dist
                        atom(i).neigh.index
                        rm_ind=[rm_ind max(atom(i).neigh.index)];
                        pause(1) %
                    elseif nNeigh < 4 % == 3
                        disp('Ti atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'Ti_un'};
                    else
                        i
                        Neighbours
                    end
                end

                % If the element is H
                if strncmp(atom(i).type,{'H'},1) % H

                    if size(Neigh_cell,1) == 1
                        atom(i).fftype={'H'};
                    elseif nNeigh > 1
                        disp('H atom over coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'H_ov'};
                        atom(i).fftype={'H'};
                        rm_ind=[rm_ind i];
                        pause(.1) %
                    elseif nNeigh < 1
                        disp('H atom under coordinated')
                        i
                        Neighbours
                        atom(i).fftype={'H_un'};
                    end
                end

                % If the element is O
                if strncmp(atom(i).type,{'O'},1)
                    if strcmp(Neighbours, 'AlAlAl') || strcmp(Neighbours, 'AlAlAlAl')
                        atom(i).fftype = {'Ob'};
                    elseif strcmp(Neighbours, 'AlAlAlAlt')
                        atom(i).fftype = {'Obt'};
                    elseif strcmp(Neighbours, 'AlAlAlH')
                        atom(i).fftype = {'Oh'};
                    elseif strcmp(Neighbours, 'AlAlAlt') || strcmp(Neighbours, 'AlAlt')
                        atom(i).fftype = {'Ops'};
                    elseif strcmp(Neighbours, 'AlAlFe') || strcmp(Neighbours, 'AlAlFeo')
                        atom(i).fftype = {'Ob'};
                    elseif strcmp(Neighbours, 'AlAlH')
                        atom(i).fftype = {'Oh'};
                    elseif strcmp(Neighbours, 'AlAlSi')
                        atom(i).fftype = {'Op'};
                    elseif strcmp(Neighbours, 'AlAlSiSi')
                        atom(i).fftype = {'Oz'};
                    elseif strcmp(Neighbours, 'AlAleH') || strcmp(Neighbours, 'AleH')
                        atom(i).fftype = {'Oh'};
                    elseif strcmp(Neighbours, 'AlAleSi') || strcmp(Neighbours, 'AleSi')
                        atom(i).fftype = {'Ob'};
                    elseif strcmp(Neighbours, 'AlAltH')
                        atom(i).fftype = {'Ops'};
                    elseif strcmp(Neighbours, 'AlFeFe') || strcmp(Neighbours, 'AlFeoFeo') || strcmp(Neighbours, 'AltFeFe') || strcmp(Neighbours, 'AltFeoFeo')
                        atom(i).fftype = {'Ops'};
                    elseif strcmp(Neighbours, 'AlFeH') || strcmp(Neighbours, 'AlFeoH')
                        atom(i).fftype = {'Oh'};
                    elseif strcmp(Neighbours, 'AlFeSi') || strcmp(Neighbours, 'AlFeoSi')
                        atom(i).fftype = {'Op'};
                    elseif strcmp(Neighbours, 'AlFet') || strcmp(Neighbours, 'AlAlFet')
                        atom(i).fftype = {'Ops'};
                    elseif strcmp(Neighbours, 'AlH')
                        atom(i).fftype = {'Oalh'}; % Al-O-H or Al-O-Si
                    elseif strcmp(Neighbours, 'AlHH')
                        atom(i).fftype = {'Oalhh'};
                    elseif strcmp(Neighbours, 'AlHSi')
                        atom(i).fftype = {'Oahs'}; % Al-OH-Si for acidic edge
                    elseif strncmp(Neighbours, 'AlHMg', 5)
                        atom(i).fftype = {'Ohmg'};
                    elseif strcmp(Neighbours, 'AlMgSi') || strcmp(Neighbours, 'AlMgoSi')
                        atom(i).fftype = {'Omg'};
                    elseif strcmp(Neighbours, 'AlOmg')
                        atom(i).fftype = {'Odsub'};
                    elseif strcmp(Neighbours, 'AltH')
                        atom(i).fftype = {'Oh'};
                    elseif strncmp(Neighbours, 'AltMgh', 6) || strcmp(Neighbours, 'AltMgoMgoMgo')
                        atom(i).fftype = {'Ops'};
                    elseif strcmp(Neighbours, 'AltSi')
                        atom(i).fftype = {'Obs'};
                    elseif strcmp(Neighbours, 'CaCaCaCaCaCa') || strcmp(Neighbours, 'CaoCaoCaoCaoCaoCao')
                        atom(i).fftype = {'Ob'};
                    elseif strcmp(Neighbours, 'CaCaCaH') || strcmp(Neighbours, 'CahCahCahH') || strcmp(Neighbours, 'CaoCaoCaoH')
                        atom(i).fftype = {'Oh'};
                    elseif strcmp(Neighbours, 'Fe2Fe2Fe2Fe2Fe2Fe2') || strcmp(Neighbours, 'FeFeFeFeFeFe')
                        atom(i).fftype = {'Op'};
                    elseif strcmp(Neighbours, 'FeFe') || strcmp(Neighbours, 'FeoFeo')
                        atom(i).fftype = {'Ob'};
                    elseif strcmp(Neighbours, 'FeFeFe') || strcmp(Neighbours, 'FeoFeoFeo')
                        atom(i).fftype = {'Ob'};
                    elseif strcmp(Neighbours, 'FeFeFet') || strcmp(Neighbours, 'FeoFeoFet')
                        atom(i).fftype = {'Ob'};
                    elseif strcmp(Neighbours, 'FeFeFeFet') || strcmp(Neighbours, 'FeoFeoFeoFet')
                        atom(i).fftype = {'Obt'};
                    elseif strcmp(Neighbours, 'FeFeFeFe') || strcmp(Neighbours, 'FeoFeoFeoFeo')
                        atom(i).fftype = {'Ob'};
                    elseif strcmp(Neighbours, 'FeFeFeH') || strcmp(Neighbours, 'FeoFeoFeoH')
                        atom(i).fftype = {'Oh'};
                    elseif strcmp(Neighbours, 'FeFeH') || strcmp(Neighbours, 'FeoFeoH')
                        atom(i).fftype = {'Oh'};
                    elseif strcmp(Neighbours, 'FeFeSi') || strcmp(Neighbours, 'FeoFeoSi')
                        atom(i).fftype = {'Op'};
                    elseif strcmp(Neighbours, 'FeH') || strcmp(Neighbours, 'FeoH')
                        atom(i).fftype = {'Oh'};
                    elseif strncmp(Neighbours, 'FeHMg', 5) || strncmp(Neighbours, 'FeoHMg', 5)
                        atom(i).fftype = {'Ohmg'};
                    elseif strcmp(Neighbours, 'FeMgSi') || strncmp(Neighbours, 'FeMgoSi', 5) || strncmp(Neighbours, 'FeoMgSi', 5)
                        atom(i).fftype = {'Omg'};
                    elseif strcmp(Neighbours, 'FeSi') || strcmp(Neighbours, 'FeoSi') || strcmp(Neighbours, 'FetSi')
                        atom(i).fftype = {'Oalt'};
                    elseif strcmp(Neighbours, 'FetFet') || strcmp(Neighbours, 'FetFetH') || strcmp(Neighbours, 'FetH')
                        atom(i).fftype = {'Oh'};
                    elseif strcmp(Neighbours, 'HLiMgMg') || strcmp(Neighbours, 'HLiMgoMgo') || strcmp(Neighbours, 'HLioMghMgh')
                        atom(i).fftype = {'Ohli'};
                    elseif strcmp(Neighbours, 'HHMg') || strcmp(Neighbours, 'HHMgh') || strncmp(Neighbours, 'HHMgo', 5)
                        atom(i).fftype = {'Omhh'};
                    elseif strcmp(Neighbours, 'HMg') || strcmp(Neighbours, 'HMgh')
                        atom(i).fftype = {'Ome'};
                    elseif strcmp(Neighbours, 'HMgMg') || strcmp(Neighbours, 'HMghMgh')
                        atom(i).fftype = {'Ohmg'};
                    elseif strcmp(Neighbours, 'HMgMgMg') || strcmp(Neighbours, 'HMghMghMgh')
                        atom(i).fftype = {'Oh'};
                    elseif strcmp(Neighbours, 'HSi')
                        atom(i).fftype = {'Osih'};
                    elseif strncmp(Neighbours, 'LiLiLiLi', 8)
                        atom(i).fftype = {'Oli'};
                    elseif strcmp(Neighbours, 'LiMgMgSi') || strcmp(Neighbours, 'LioMgMgSi') || strcmp(Neighbours, 'LioMgoMgoSi')
                        atom(i).fftype = {'Oli'};
                    elseif strcmp(Neighbours, 'MgMgMgMgMgMg') || strcmp(Neighbours, 'MghMghMghMghMghMgh') || strcmp(Neighbours, 'MgoMgoMgoMgoMgoMgo')
                        atom(i).fftype = {'Ob'};
                    elseif strcmp(Neighbours, 'MgMgMgSi') || strcmp(Neighbours, 'MghMghMghSi') || strcmp(Neighbours, 'MgoMgoMgoSi')
                        atom(i).fftype = {'Op'};
                    elseif strcmp(Neighbours, 'MgMgSi') || strcmp(Neighbours, 'MgoMgoSi')
                        atom(i).fftype = {'Odsub'};
                    elseif strcmp(Neighbours, 'MgSi')
                        atom(i).fftype = {'Omg'};
                    elseif strcmp(Neighbours, 'SiSi') || strcmp(Neighbours, 'SiSiSi')
                        atom(i).fftype = {'Ob'};
                    elseif strcmp(Neighbours, 'TiTiTi') || strcmp(Neighbours, 'TioTioTio')
                        atom(i).fftype = {'Ob'};
                        % Handle special case for 'AlSi' with nested conditions
                    elseif strcmp(Neighbours, 'AlSi')
                        if sum(ismember(find(strcmp([atom.fftype], 'Alt')), [atom(i).neigh.index])) > 0
                            atom(i).fftype = {'Oalt'};
                        else
                            atom(i).fftype = {'Oalsi'}; % Special zeolite case
                        end %% END on new statement
                    elseif strcmp(Neighbours,'HH')
                        atom(i).fftype={'Ow'};
                        disp('Water?')
                    elseif nNeigh > 2
                        disp('O atom overcoordinated')
                        i
                        Neighbours
                        atom(i).fftype={'O_ov'};
                    elseif nNeigh == 1 || strcmp(Neighbours,'AlAl') || strcmp(Neighbours,'AlMg') || strcmp(Neighbours,'Si')
                        if strcmp(Neighbours,'Si')
                            atom(i).fftype={'Osi'};
                        elseif strncmp(Neighbours,'AlAl',2)
                            atom(i).fftype={'Oal'};
                        else
                            disp('O atom undercoordinated')
                            i
                            Neighbours
                            atom(i).fftype={'O_un'};
                        end
                    elseif nNeigh == 0
                        Neighbours
                        atom(i)
                        atom(i).neigh.dist
                        %                     atom(i)=[];
                        disp('remove O atom...?')

                        rm_ind=[rm_ind i];
                        pause(1) %
                    end

                    if strcmp(atom(i).fftype,'Oalh') && Add_H == 1 % && assignment_run == 8
                        disp('Adding acidic H to Oalh-->Oalhh')
                        atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
                        atom(end).index=size(atom,2);
                        Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.25);
                        NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) -1/3*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1/3*Neigh.r_vec(:,3)],1)));
                        [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                    end
                end
                All_Neighbours=[All_Neighbours;{char(atom(i).type) char(atom(i).fftype) Neighbours}]; % New 2.11
            end
        i=i+1;
        [atom.type]=atom.fftype;
    end
end

Atom_labels=unique([atom.element]);
temp=element_atom(atom);
[atom.element]=temp.type;

All_Neighbours=sortrows(All_Neighbours,1);
All_Neighbours=sortrows(All_Neighbours,3);
All_Neighbours=sortrows(All_Neighbours,1);

i=1;All_Neighbours(:,4)={1};
while i<size(All_Neighbours,1)+1
    if i==1
        All_Neighbours{i,4}=1;
    elseif strcmp(All_Neighbours(i,1),All_Neighbours(i-1,1)) && strcmp(All_Neighbours(i,3),All_Neighbours(i-1,3))
        All_Neighbours{i-1,4}=All_Neighbours{i-1,4}+1;
        All_Neighbours(i,:)=[];
        i=i-1;
    else
        All_Neighbours{i,4}=1;
    end
    i=i+1;
end

All_Neighbours=[All_Neighbours(:,1) All_Neighbours];
i=1;
while i<size(All_Neighbours,1)+1
    n=0;
    if sum(strcmp(All_Neighbours(i,1),All_Neighbours(:,1))) > 1
        n=1;
    else
        n=0;
    end

    if i>1 && sum(strcmp(All_Neighbours(i,1),All_Neighbours(:,1))) > 1
        n=sum(strcmp(All_Neighbours(i,1),All_Neighbours(1:i,1)));
    end

    if n>0
        All_Neighbours{i,2}=strcat(All_Neighbours{i,2},num2str(n));
    else

end
All_Neighbours(i,3)
    All_Neighbours(i,6)={unique([atom(strcmp([atom.type],All_Neighbours(i,3))).charge])};
    i=i+1;
end

assignin('caller','remove_ind',rm_ind);
assignin('caller','All_Neighbours',All_Neighbours);


