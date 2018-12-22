%% find_bonded_atom.m
% * This function tries to find all bonds between the  atomtype1 and
% * atomtype2.
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # fb_atom = find_bonded_atom(atom,Box_dim,{'Oalhh'})
% # fb_atom = find_bonded_atom(atom,Box_dim,)[1:4:100])
% # fb_atom = find_bonded_atom(atom,Box_dim,atomtype1,atomtype2)
% # fb_atom = find_bonded_atom(atom,Box_dim,atomtype1,atomtype2,1.25,2.25)

function fb_atom = find_bonded_atom(atom,Box_dim,atomtype1,varargin)

if isnumeric(atomtype1)
    indtype1=atomtype1;
    atomtype1=[atom(indtype1(1)).type];
else
    indtype1=find(ismember([atom.type],atomtype1))';
end

if nargin > 3
    atomtype2=varargin{1};
    indtype2=find(ismember([atom.type],atomtype2))';
else
    indtype2=[1:size(atom,2)]';
    indtype2=setdiff(indtype2,indtype1);
    atomtype2=unique([atom(indtype2).type]);
end


if nargin > 4
    rmaxshort=varargin{2};
    rmaxlong=varargin{3};
else
    rmaxshort=1.25;
    rmaxlong=2.25;
end

if numel(indtype1)==0
    try
        indtype1=find(strcmp([atom.type],atomtype1));
    catch
        indtype1=find(strncmpi([atom.type],atomtype1,1));
    end
end

if numel(indtype2)==0
    try
        indtype2=find(strcmp([atom.type],atomtype2));
    catch
        try
        indtype2=find(strncmpi([atom.type],atomtype2,1));
        catch
            disp('Matlab cannot fetch the indexes for all the atomtypes given.')
        end
    end
end

%dist_matrix = cell_list_dist_matrix_atom(atom,Box_dim,rmaxshort,rmaxlong);
dist_matrix = dist_matrix_atom(atom,Box_dim);

Bond_index_long=[];
indtype12long=[];
for i=1:numel(indtype1)
    neigh_index=find(dist_matrix(:,indtype1(i))<rmaxlong&dist_matrix(:,indtype1(i))>0);
    neigh_index=intersect(neigh_index,indtype2);
    indtype12long=[indtype12long;neigh_index];
    Bond_index_long=[Bond_index_long; [indtype1(i).*ones(numel(neigh_index),1) neigh_index dist_matrix(neigh_index,indtype1(i))]];
end
indtype12long=sort(indtype12long);
indtype12long=setdiff(indtype12long,indtype1)';
if size(indtype12long,2)>size(indtype12long,1)
    indtype12long=indtype12long';
end

Bond_index_short=[];
indtype12short=[];
for i=1:numel(indtype1)
    neigh_index=find(dist_matrix(:,indtype1(i))<rmaxshort&dist_matrix(:,indtype1(i))>0);
    neigh_index=intersect(neigh_index,indtype2);
    indtype12short=[indtype12short;neigh_index];
    Bond_index_short=[Bond_index_short; [indtype1(i).*ones(numel(neigh_index),1) neigh_index dist_matrix(neigh_index,indtype1(i))]];
end
% indtype12short=sort(indtype12short);
% indtype12short=setdiff(indtype12short,indtype1)';
% if size(indtype12short,2)>size(indtype12short,1)
%     indtype12short=indtype12short';
% end

if strncmpi(atomtype1,'H',1)
%     indtype2=indtype12short;
    Bond_index=Bond_index_short;
elseif strncmpi(atomtype2,'H',1)
%     indtype2=indtype12short;
    Bond_index=Bond_index_short;
else
%       indtype2=[indtype12long; indtype12short];
    Bond_index=[Bond_index_long;Bond_index_short];
end

[Y,I]=sort(Bond_index(:,1));
Bond_index=Bond_index(I,:);
Bond_index = unique(Bond_index,'rows','stable');

indtype1=unique(Bond_index(:,1));
indtype2=unique(Bond_index(:,2));
if size(indtype1,2)>size(indtype1,1)
    indtype1=indtype1';
end
if size(indtype2,2)>size(indtype2,1)
    indtype2=indtype2';
end

fb_atom=atom(sort([indtype1; indtype2]));

assignin('caller','dist_matrix',dist_matrix);
assignin('caller','Bond_index',Bond_index);
assignin('caller','type1_atom',atom(sort(indtype1)));
assignin('caller','type2_atom',atom(sort(indtype2)));
assignin('caller','type1_ind',sort(indtype1));
assignin('caller','type2_ind',sort(indtype2));

