%% ff_atom.m
% * This function tries to assign all atoms according to some custom force field
%
%% Version
% 3.0.2
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom=ff_atom(atom,Box_dim)
% # atom=ff_atom(atom,Box_dim,rmaxlong)
% # atom=ff_atom(atom,Box_dim,rmaxlong,Bond_index)


%% See also
% # <adjust_H_atom.html adjust_H_atom>
% # <charge_minff_atom.html charge_minff_atom>
% # <cn_atom.html cn_atom>
% # <element_atom.html element_atom>
% # <number_type.html number_type>
%
function atom=ff_atom(atom,Box_dim,varargin)
%%
format compact;

rmaxlong=2.45;
if nargin>2
    rmaxlong=varargin{1};
end

if nargin>3
    Bond_index=varargin{2};
end

temp_atom=atom;

% % atom=element_atom(atom);
% % [atom.element]=atom.type;
% % [atom.fftype]=atom.element;
% %
% % XYZ_labels=[atom.type]';
% % XYZ_data=[[atom.x]' [atom.y]' [atom.z]'];

%% From minff_atom
[atom.element]=atom.type;

for i=1:size(atom,2)
    if strncmpi(atom(i).type,{'Si'},2);atom(i).element={'Si'};
    elseif strncmp(atom(i).type,{'SC'},2);atom(i).element={'Si'};
    elseif strncmpi(atom(i).type,{'Ale'},3);atom(i).element={'Ale'};
    elseif strncmpi(atom(i).type,{'Alt'},3);atom(i).element={'Alt'};
    elseif strncmpi(atom(i).type,{'Al'},2);atom(i).element={'Al'};
    elseif strncmpi(atom(i).type,{'Mge'},3);atom(i).element={'Mge'};
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
%% From minff_atom

if exist('Bond_index','var')
    atom=cn_atom(atom,Box_dim,rmaxlong,Bond_index);
else
    atom=cn_atom(atom,Box_dim,rmaxlong);
end

nBonds=size(Bond_index,1);
assignin('caller','nBonds',nBonds);

assignin('caller','Bond_index',Bond_index);
assignin('caller','dist_matrix',dist_matrix);

try
    assignin('caller','Neigh_index',Neigh_index);
    assignin('caller','radius_limit',radius_limit);
catch
    disp('Could not assign Neigh_index or radius_limit')
end

All_Neighbours=[];All_types=struct;
Heal_O=0;
i=1;
while i <= size(atom,2)
    if mod(i,100)==1
        i-1
    end

    if strncmpi([atom(i).resname],'SOL',3)==0 && strncmpi([atom(i).resname],'ION',3)==0
        nNeigh=numel([atom(i).neigh]);
        if nNeigh>0
            Neigh_cell = sort([atom(i).neigh.type]);
            if length(Neigh_cell) > 0
                Neighbours=strcat(Neigh_cell{:});
                Neighboursx1=cellfun(@(x) x(1),Neigh_cell,'un',0);
                Neighboursx1=lower(strcat(Neighboursx1{:}));
                All_Neighbours=[All_Neighbours;{char(atom(i).type) Neighbours Neighboursx1 } i]; %
            else
                Neighbours={'Nan'};
                Neighboursx1={'Nan'};
            end
        end
    end
    i=i+1;
    % [atom.type]=atom.fftype;
end

if Heal_O>0
    atom=adjust_H_atom(atom,Box_dim);
end

if numel(All_Neighbours)>0

    All_Neighbours=sortrows(All_Neighbours,1);
    All_Neighbours=sortrows(All_Neighbours,2);
    All_Neighbours=sortrows(All_Neighbours,1);

    i=1;
    % All_Neighbours(:,5)=All_Neighbours(:,4);
    All_Neighbours(:,5)={1};
    while i<size(All_Neighbours,1)+1
        if i==1
            All_Neighbours{i,5}=1;
        elseif strcmp(All_Neighbours(i,1),All_Neighbours(i-1,1)) && strcmp(All_Neighbours(i,2),All_Neighbours(i-1,2))
            All_Neighbours{i-1,5}=All_Neighbours{i-1,5}+1;
            All_Neighbours{i-1,4}=[All_Neighbours{i-1,4} All_Neighbours{i,4}];
            All_Neighbours(i,:)=[];
            i=i-1;
        else
            All_Neighbours{i,5}=1;
        end
        i=i+1;
    end

    %% Set fftype names with numbers
    All_Neighbours=[All_Neighbours All_Neighbours(:,1)];
    % All_Neighbours(:,2)=lower(All_Neighbours(:,2));
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
            All_Neighbours{i,6}=strcat(All_Neighbours{i,6},num2str(n));
        else

        end

        i=i+1;
    end

    All_Neighbours=[All_Neighbours All_Neighbours(:,1) ];
    % All_Neighbours(:,1)=regexprep(All_Neighbours(:,1),'\d+$','');

    for i=1:size(All_Neighbours,1)
        ind=All_Neighbours{i,4};
        if ~ismember(All_Neighbours{i,1},{'OW' 'Ow' 'HW1' 'HW2' 'HW' 'Hw'})
            if ismember(All_Neighbours{i,1},{'O' 'H'})
                for j=1:size(ind,2)
                    % if strncmp(All_Neighbours{i,1},'H',1)
                    %     [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),All_Neighbours{i,3})
                    % elseif length(All_Neighbours{i,3})==1
                    %     [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'s',All_Neighbours{i,3});
                    % elseif length(All_Neighbours{i,3})==2
                    %     [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'b',All_Neighbours{i,3});
                    % elseif length(All_Neighbours{i,3})==3
                    %     [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'p',All_Neighbours{i,3});
                    % elseif length(All_Neighbours{i,3})==4
                    %     [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'p',All_Neighbours{i,3});
                    % elseif length(All_Neighbours{i,3})==5
                    %     [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'e',All_Neighbours{i,3});
                    % elseif length(All_Neighbours{i,3})==6
                    %     [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'o',All_Neighbours{i,3});
                    % end
                    % % else
                    % %     [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),All_Neighbours(i,3));
                    if strncmp(All_Neighbours{i,1},'H',1)
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),All_Neighbours{i,3})
                    elseif length(All_Neighbours{i,3})==1
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),All_Neighbours{i,3});
                    elseif length(All_Neighbours{i,3})==2
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),All_Neighbours{i,3});
                    elseif length(All_Neighbours{i,3})==3
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),All_Neighbours{i,3});
                    elseif length(All_Neighbours{i,3})==4
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),All_Neighbours{i,3});
                    elseif length(All_Neighbours{i,3})==5
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),All_Neighbours{i,3});
                    elseif length(All_Neighbours{i,3})==6
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),All_Neighbours{i,3});
                    end
                end
            else
                for j=1:size(ind,2)
                    if length(All_Neighbours{i,3})==1
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'d');
                    elseif length(All_Neighbours{i,3})==2
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'b');
                    elseif length(All_Neighbours{i,3})==3
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'p');
                    elseif length(All_Neighbours{i,3})==4
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'t');
                    elseif length(All_Neighbours{i,3})==5
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'e');
                    elseif length(All_Neighbours{i,3})==6
                        [atom(ind(j)).fftype]=strcat(All_Neighbours(i,1),'o');
                    else
                        [atom(ind(j)).fftype]=All_Neighbours(i,1);
                    end
                end
            end
        end
    end

    All_Neighbours=[All_Neighbours All_Neighbours(:,1)];
    for i=1:size(All_Neighbours,1)
        ind=All_Neighbours{i,4};
        All_Neighbours{i,7}=temp_atom(ind(1)).type;
        All_Neighbours{i,8}=atom(ind(1)).fftype;
    end


    if ~isfield(atom,'charge')
        atom=charge_minff_atom(atom,Box_dim,{'Al' 'Alo' 'Alt' 'Ale' 'Tio' 'Feo3' 'Fet3' 'Fee3' 'Feo2' 'Fet2' 'Fee2' 'Fs' 'Na' 'K' 'Cs' 'Mgo' 'Mgh' 'Mge' 'Cao' 'Cah' 'Ca' 'Cl' 'Sit' 'Si' 'Sio' 'Site' 'Lio' 'H'},...
            [1.782 1.782 1.782 1.985 2.48 1.5 1.5 1.75 1.184 1.184 1.32 -0.76 1 1 1 1.562 1.74 1.635 1.66 1.52 2 -1 1.884 1.884 1.884 2.413 0.86 0.4]);
    end

    if isfield(atom,'charge')
        for i=1:size(All_Neighbours,1)
            ind=All_Neighbours{i,4};
            All_Neighbours(i,9)={unique(round2dec([atom(ind).charge],8))};
        end
    end

    if isfield(atom,'cn')
        for i=1:size(All_Neighbours,1)
            ind=All_Neighbours{i,4};
            All_Neighbours(i,10)={unique([atom(ind).cn])};
        end
    end

    % atom = number_type(atom,0);
    % ff=load(strcat(ffname,'_ff.mat'));
    % ff=ff.ff;
    % [ff.resname]=deal({'MIN'});
    % ff=resname_atom(ff);
    % ff(ismember([ff.resname],{'ION' 'SOL'}))=[];
    %
    % for i=1:size(ff,2)
    %     [ff(i).type]=lower([ff(i).type]);
    % end
    % ff=ff(ismember([ff.type],unique(All_Neighbours(:,1))));

    % % i=1;
    % % while i<size(ff,2)+1
    % %     ind=find(strncmpi([atom.type],[ff(i).type],2));
    % %     if numel(ind)==0
    % %         ind=find(strncmpi([atom.type],[ff(i).type],1));
    % %     end
    % %     if numel(ind)>0
    % %         for j=1:size(ind,2)
    % %             [atom(ind(j)).charge]=[ff(i).charge];
    % %             [atom(ind(j)).e_kJmol]=[ff(i).e_kJmol];
    % %             [atom(ind(j)).sigma_nm]=[ff(i).sigma_nm];
    % %         end
    % %         i=i+1;
    % %     else
    % %         ff(i)=[];
    % %     end
    % % end
    %
    % assignin('caller','ff',ff);

end

[atom.type]=atom.fftype;

for i=1:length(unique([atom.type]))
    new_Atom_label=sort(unique([atom.type]));
    try
        ind=ismember([atom.type],new_Atom_label(i));
        assignin('caller',strcat(char(new_Atom_label(i)),'_atom'),atom(ind));
    catch
        disp('Could not finalize:')
        i
        new_Atom_label
    end
end

for i=1:size(All_Neighbours,1)
    All_types(i).type=All_Neighbours{i,7};
    All_types(i).typenum=All_Neighbours(i,6);
    All_types(i).element=All_Neighbours(i,1);
    All_types(i).neighbours=All_Neighbours(i,2);
    All_types(i).neigh=All_Neighbours(i,3);
    All_types(i).index=All_Neighbours{i,4};
    All_types(i).number=All_Neighbours{i,5};
    All_types(i).type2=All_Neighbours{i,8};
    All_types(i).charge=All_Neighbours{i,9};
    All_types(i).cn=All_Neighbours{i,10};
end

assignin('caller','Bond_index',Bond_index);
assignin('caller','All_Neighbours',All_Neighbours);
assignin('caller','All_types',All_types);
assignin('caller','XYZ_labels',XYZ_labels);
assignin('caller','XYZ_data',XYZ_data);


