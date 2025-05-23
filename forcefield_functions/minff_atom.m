%% minff_atom.m
% * This function tries to assign all atoms according to the MINFF atom
% types.
% * heal_iterations should be a list of numbers 1-7, corresponding to which
% assigment runs you need in order to heal Al, Mg, Si, O, H and so on...
% * The variables distance_factor and rmaxlong are related to the
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
% # atom=minff_atom(atom,Box_dim,'minff')
% # atom=minff_atom(atom,Box_dim,'minff',heal_iterations)

function atom=minff_atom(atom,Box_dim,varargin)
%%
format compact;

if nargin>2
    ffname=varargin{1};
else
    ffname='minff';
end

if nargin>3
    heal_iterations=varargin{2};
else
    heal_iterations=[1 2]; % No healing, will run much faster
end

distance_factor=0.6;
rmaxlong=2.45;

if nargin>4
    distance_factor=varargin{3};
end
if nargin>5
    rmaxlong=varargin{4};
end

atom=element_atom(atom);

% Initialize some variables
rm_ind=[];
Heal_Al=0;
Heal_Mgo=0;
Heal_Feo=0;
Heal_Si=0;
Heal_O=0;
Heal_H=0; % I do not think this is important
Add_H=0; % Add H to get neutral edge
Add_extra_H=0; % Add H to get positive edge As a final step, add extra H's to a Oalhh atoms...nH_extra=16;
% n=7; % with edge structure
for assignment_run=heal_iterations

    All_Neighbours=[];
    % Heal broken bonds in the structure in steps, start from the center of the particle
    if assignment_run==1
        Heal_Al=0; % 1 for yes and 0 for no
    elseif assignment_run==2
        Heal_Mgo=0; % 1 for yes and 0 for no
    elseif assignment_run==3
        Heal_Feo=0; % 1 for yes and 0 for no
    elseif assignment_run==4
        Heal_Si=1; % 1 for yes and 0 for no
    elseif assignment_run==5
        Heal_O=1; % 1 for yes and 0 for no
    elseif assignment_run==6
        Heal_H=0; % 1 for yes and 0 for no
    elseif assignment_run==7
        Add_H=1 % 1 for yes and 0 for no
    elseif assignment_run==8
        Add_extra_H=1 % As a final step, add extra H's to a Oalhh atoms...nH_extra=16;
    end

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
                        if Heal_Si == 1
                            atom(size(atom,2)+1)=atom(find(strcmp([atom.type],'O'),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+1.61*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1.5*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1.5*Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
                    else
                        i
                        Neighbours
                    end
                end

                % If the element is Al
                if strncmpi(atom(i).type,{'Al'},2) % Al

                    if sum(strncmp({'O'},[atom(i).neigh.type],1)) == 6 % Al O O O O O O
                        atom(i).fftype={'Alo'};
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
                        if Heal_Al == 1
                            disp('Healing Al with O!!!')
                            atom(size(atom,2)+1)=atom(find(strncmpi([atom.type],'O',1),1));
                            atom(end).index=size(atom,2);
                            atom(end).type={'O'};
                            atom(end).fftype={'O'};
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]-1.85*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
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

                        if mean([atom(i).neigh.dist])>2.05 && mean([atom(i).neigh.dist]<2.1)
                            disp('Ave oct Fe-O distance is large..')
                            mean([atom(i).neigh.dist])
                            if mean([atom(i).neigh.dist])<2.07
                                disp('Guessing its Fe3+..')
                            else
                                disp('Guessing its Fe2+..')
                            end
                        end

                        if mean([atom(i).neigh.dist])<2.07
                            atom(i).fftype={'Feo3'};
                        else
                            atom(i).fftype={'Feo2'};
                        end
                        
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 5 % Fe O O O O
                        atom(i).fftype={'Fee3'};
                    elseif sum(strncmp({'O'},[atom(i).neigh.type],1)) == 4 % Fe O O O O
                        if sum([atom(i).neigh.dist])<7.5
                            atom(i).fftype={'Fet3'};
                        else
                            atom(i).fftype={'Fet2'};
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
                        if Heal_Feo == 1
                            atom(size(atom,2)+1)=atom(find(strncmpi([atom.type],'O',1),1));
                            atom(end).index=size(atom,2);
                            atom(end).type={'O'};
                            atom(end).fftype={'O'};
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]-1.9*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        end
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
                        if Heal_H == 1
                            atom(size(atom,2)+1)=atom(find(strncmpi([atom.type],'O',1),1));
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.6);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y atom(i).z]+mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
                        elseif nNeigh < 1
                            Neighbours
                            %                             atom(i)=[];
                            %                             disp('removed atom...')
                            i
                            pause(0.1) %
                        end
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
                        if Heal_O == 1 % && assignment_run == 6
                            disp('Adding H to unsaturated O')
                            try
                                atom(size(atom,2)+1)=atom(find(strncmp([atom.type],'H',1),1));
                            catch
                                atom(size(atom,2)+1)=atom(i);
                            end
                            atom(end).type={'H'};
                            atom(end).fftype={'H'};
                            atom(end).index=size(atom,2);
                            Neigh = neighbor_func(i,[[atom(i).x]' [atom(i).y]' [atom(i).z]'],[[atom.x]' [atom.y]' [atom.z]'],Box_dim,2.25);
                            NewNeighCoords=num2cell([atom(i).x atom(i).y (atom(i).z)]-1*mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1*Neigh.r_vec(:,3)],1)/norm(mean([Neigh.r_vec(:,1) Neigh.r_vec(:,2) 1*Neigh.r_vec(:,3)],1)));
                            [atom(end).x atom(end).y atom(end).z]=deal(NewNeighCoords{:});
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
                All_Neighbours=[All_Neighbours;{char(atom(i).type) char(atom(i).fftype) Neighbours}]; % New 2.11
            end
        end
        i=i+1;
        [atom.type]=atom.fftype;
    end
end

if assignment_run>7
    atom=adjust_H_atom(atom,Box_dim);
end

Atom_labels=unique([atom.element]);
temp=element_atom(atom);
[atom.element]=temp.type;

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

atom=charge_minff_atom(atom,Box_dim,{'Al' 'Alo' 'Alt' 'Ale' 'Tio' 'Feo3' 'Fet3' 'Fee3' 'Feo2' 'Fet2' 'Fee2' 'Fs' 'Na' 'K' 'Cs' 'Mgo' 'Mgh' 'Mge' 'Cao' 'Cah' 'Sit' 'Si' 'Sio' 'Site' 'Lio' 'H'},[1.782 1.782 1.782 1.985 2.48 1.5 1.5 1.75 1.184 1.184 1.32 -0.76 1 1 1 1.562 1.74 1.635 1.66 1.52 1.884 1.884 1.884 2.413 0.86 0.4]);

if abs(round2dec(sum([atom.charge])) - sum([atom.charge]))>0.0001
    disp('Initial total charge!')
    Total_charge=round2dec(sum([atom.charge]),8)
    [atom(1).charge]=atom(1).charge+(round2dec(sum([atom.charge])) - sum([atom.charge]));
    %     atom=tweak_charge_atom(atom,'Al');
    disp('Final total charge, modfying the first atoms charge somewhat!')
    Total_charge=sum([atom.charge])
    round2dec(sum(Total_charge),6)
end

for i=1:length(unique([atom.type]))
    try
        new_Atom_label=sort(unique([atom.type]));
        ind=ismember([atom.type],new_Atom_label(i));
        assignin('caller',strcat(char(new_Atom_label(i)),'_atom'),atom(ind));
    catch
        disp('Could not finalize:')
        i
        Atom_labels
    end
end

atom = order_attributes(atom);

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

% atom_ff=ff_atom(atom,Box_dim,'minff',rmaxlong);
% assignin('caller','atom_ff',atom_ff);
% assignin('caller','All_types',All_types);

assignin('caller','remove_ind',rm_ind);
assignin('caller','All_Neighbours',All_Neighbours);


