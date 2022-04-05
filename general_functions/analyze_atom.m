%% analyze_atom.m
% * This function fetches the ionic radius, originally taken from below
% * Ref. 1	Revised effective ionic radii and systematic studies of
% interatomic distances in halides and chalcogenides. R. D. Shannon Acta
% Cryst. (1976) A32, 751-767.
% * Ref. 2	Electronic Table of Shannon Ionic Radii, J. David Van Horn,
% 2001, downloaded MO/DA/YEAR.
% *
% * This function also calculates the bond valence values according to
% http://www.iucr.org/resources/data/datasets/bond-valence-parameters
% compiled by I. David Brown, McMaster University, Ontario, Canada
% * Data set bvparm2016.cif: 2016 version, (posted 2016-11-03)
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # analyze = analyze_atom(atom,Box_dim)
% # analyze = analyze_atom(atom,Box_dim,2.5)
% # analyze = analyze_atom(atom,Box_dim,2.5,Bond_index,Valencestates)
%
function properties = analyze_atom(atom,Box_dim,varargin)

% FIX multiple OxStates and other from the Shannon file
load('Revised_Shannon_radii.mat');

if nargin > 2
    rmax=varargin{1};
else
    rmax=2.25;
end

element=element_atom(atom); % Optionally, set the names ow water O and water H
% element=bond_angle_atom(element,Box_dim,1.25,rmax,'more');

if nargin < 4
    disp('Trying to find bonded atoms')
    element=bond_atom(element,Box_dim,rmax);
else
    disp('Trying to recalculate the bonded atoms')
    Bond_index=varargin{2};
    element=recalc_bond_atom(element,Box_dim,Bond_index);
end

disp('Running Bond Valence analysis...')
if nargin > 4
    Valencestates=varargin{3};
    element=bond_valence_atom(element,Box_dim,1.25,rmax,Valencestates);
else
    element=bond_valence_atom(element,Box_dim,1.25,rmax);
end

Valences=[element.valence];
modValences=sum(Valences-round(Valences))/numel(Valences);
% other_ind=[];ValenceGuess=zeros(1,size(element,2));
% for i=1:size(element,2)
%     try
%         element(i).type;
%         ValenceGuess(i)=round([element(i).valence]);
%         element(i)=bond_valence_atom(element(i),Box_dim,1.25,rmax,ValenceGuess(i));
%     catch
%         other_ind=[other_ind i];
%     end
%     if mod(i,100)==1
%         i-1
%     end
% end
% if numel(other_ind)>0
%     element(other_ind)=bond_valence_atom(element(other_ind),Box_dim,1.25,rmax);
% end
% assignin('caller','ValenceGuess',ValenceGuess)
% assignin('caller','other_ind',other_ind)

element=mass_atom(element,Box_dim);
assignin('caller','element',element);
assignin('caller','Box_volume',Box_volume);
assignin('caller','Box_density',Box_density);

for i=1:size(element,2)
    Ion_ind=find(strcmp([element(i).type],Ion));
    if numel(Ion_ind)==0
        Ion_ind=find(strncmpi([element(i).type],Ion,1));
    end
    CN_ind=find(numel(element(i).neigh.index)==CN);
    Ox_ind=find(round(element(i).valence-modValences)==OxState);
    CN_ind=intersect(Ion_ind,CN_ind);
    Ox_ind=intersect(Ion_ind,Ox_ind);
    ind=intersect(Ox_ind,CN_ind);
    
    if numel(ind)==0
        if numel(Ox_ind)==0 && numel(CN_ind)>0
            ind=CN_ind(1);
        elseif numel(CN_ind)==0 && numel(Ox_ind)>0
            ind=Ox_ind(1);
        else
            ind=Ion_ind(1);
        end
    elseif numel(ind)>1
        ind=ind(1);
    end
    properties(i).index=element(i).index;
    properties(i).type=Ion(ind);
    properties(i).fftype=atom(i).type;
    properties(i).neigh=element(i).neigh;
    properties(i).bond=element(i).bond;
    %     properties(i).angle=element(i).angle; % Because bond_atom() and not bond_angle_atom()
    properties(i).bv=element(i).bv;
    properties(i).valence=element(i).valence-modValences;
    properties(i).ave_dist=mean(properties(i).bond.dist);
    properties(i).std_dist=std(properties(i).bond.dist);
    %     properties(i).exp_dist=properties(i).ave_dist+element(i).Rdiff;
    properties(i).rdiffvalence=element(i).Rdiff;
    properties(i).cn_bv=size(properties(i).bv,2);
    properties(i).ShannonParam={'>>>>'};
    properties(i).atnum=element(i).atnum;
    properties(i).mass=element(i).mass;
    properties(i).oxstate=OxState(ind);
    properties(i).cn_guess=CN(ind);
    properties(i).ionicradii=IonicRadii(ind);
    properties(i).RevShannonind=ind;
    properties(i).type=Ion(properties(i).RevShannonind);
    properties(i).MODRevShannonind=properties(i).RevShannonind;
    %         properties(i).oxstate=OxState(properties(i).RevShannonind);
    properties(i).ip=ZoverIR(properties(i).RevShannonind);
    properties(i).cn_guess=CN(properties(i).RevShannonind);
    properties(i).crysradii=CrysRadii(properties(i).RevShannonind);
    properties(i).ionicradii=IonicRadii(properties(i).RevShannonind);
    properties(i).vdwradii=radius_vdw([element(i).type]);
    properties(i).ip=ZoverIR(properties(i).RevShannonind);
    properties(i).elecconf=strtrim(ElecConf(properties(i).RevShannonind));
    properties(i).spinstate=strtrim(SpinState(properties(i).RevShannonind));
    
    if mod(i,100)==1
        i-1
    end
end

pre_properties=properties;
assignin('caller','pre_properties',pre_properties);

for i=1:size(properties,2)
    Ion_ind=find(strcmp([element(i).type],Ion));
    if numel(Ion_ind)==0
        Ion_ind=find(strncmpi([element(i).type],Ion,1));
    end
    CN_ind=find(numel(element(i).neigh.index)==CN);
    ind=intersect(Ion_ind,CN_ind);
    
    if numel(ind)>0
        current_radii=[properties(i).ionicradii];
        neigh_radii=[properties([properties(i).neigh.index]).ionicradii];
        sum_radii=repmat(current_radii',numel(neigh_radii),1)+repmat(neigh_radii,numel(current_radii),1)';
        [minvalue,preox_ind]=min(abs(mean(sum_radii-median([properties(i).neigh.dist] - properties(i).rdiffvalence))));
        properties(i).type=Ion(ind(preox_ind));
        %         properties(i).oxstate=OxState(ind(preox_ind));
        properties(i).valence=element(i).valence;
        properties(i).ip=ZoverIR(ind(preox_ind));
        properties(i).cn_guess=CN(ind(preox_ind));
        properties(i).crysradii=CrysRadii(ind(preox_ind));
        properties(i).ionicradii=IonicRadii(ind(preox_ind));
        properties(i).vdwradii=radius_vdw([element(i).type]);
        properties(i).ip=ZoverIR(ind(preox_ind)); %preox_ind);
        properties(i).elecconf=strtrim(ElecConf(ind(preox_ind)));
        properties(i).spinstate=strtrim(SpinState(ind(preox_ind)));
        properties(i).MODRevShannonind=ind(preox_ind);
    end
    if mod(i,100)==1
        i-1
    end
end

% % From Vesta manual... Test with pyrophyllite...
% % for i=1:size(properties,2)
% %     properties(i).DistBondECoN={'>>>>'};
% %     dist=[properties(i).neigh.dist];
% %     properties(i).DistorIndex=sum((dist-mean(dist))./mean(dist))/numel(dist);
% %     properties(i).lmin=min(dist);
% %     properties(i).la=sum(dist.*exp(1-(dist./properties(i).lmin).^6))/sum((exp(1-(dist./properties(i).lmin).^6)));
% %     properties(i).wi=exp(1-(dist./properties(i).la).^6);
% %     properties(i).ECoN=sum([properties(i).wi]);
% % end
% % for i=1:size(properties,2)
% %     ind=[properties(i).neigh.index];
% %     properties(i).deltaq=-([properties(ind).oxstate].*[properties(i).wi]')./[properties(ind).ECoN];
% %     properties(i).Q=sum([properties(i).deltaq]);
% % end

diff_ind=find([properties.RevShannonind]-[properties.MODRevShannonind]);
if numel(diff_ind)>0
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Possible problems with these sites due to non-ideal bonding distances,')
    disp('Compare with output in the pre_properties struct variable')
    unique([atom(diff_ind).type])
    diff_ind
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

diff_valence=(abs([properties.oxstate])-[properties.valence])';
ind=find(abs(diff_valence)>0.5);
if numel(ind)>0
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Possible problems with these sites due to multiple oxidation states,')
    disp('or atoms that are undersaturated!')
    unique([atom(ind).type])
    ind
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end
assignin('caller','heal_ind',ind');

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Total valence from oxstate');
Tot_valence_oxstate=sum([properties.oxstate])

ind_neg=find([properties.oxstate]<0);
Tot_valence=[properties.valence];
Tot_valence(ind_neg)=-Tot_valence(ind_neg);
disp('Total valence');
Tot_valence=sum(Tot_valence)

disp('Global instability index is:')
GII=(sum((abs([properties.oxstate])-[properties.valence]).^2)/size(properties,2))^0.5
if GII>0.2
    disp('GII > 0.2 --> Structure likely not super stable...');
end

GII_noH=0;
if sum(strncmpi([properties.type],'H',1))>0
    disp('Global instability (ignoring the H...) index is:')
    ind_noH=find(~strncmpi([properties.type],'H',1));
    GII_noH=(sum((abs([properties(ind_noH).oxstate])-[properties(ind_noH).valence]).^2)/numel(ind_noH))
    if GII_noH>0.2
        disp('GII_noH > 0.2 --> Structure likely not super stable...');
    end
end

Atom_labels=unique([atom.fftype]);
for i=1:length(Atom_labels)
    ind=find(strcmp([atom.fftype],Atom_labels(i)));
    BondSummary(i).type=Atom_labels(i);
    BondSummary(i).GII=(sum((abs([properties(ind).oxstate])-[properties(ind).valence]).^2)/numel(ind))^0.5;
    BondSummary(i).d_strain=sum([properties(ind).valence]-abs([properties(ind).oxstate]))/numel(ind);
    BondSummary(i).ValenceAve=mean([properties(ind).valence]);
    BondSummary(i).ValenceStd=std([properties(ind).valence]);
    BondSummary(i).dist=mean([properties(ind).ave_dist]);
end

load('bond_valence_values.mat');
for i=1:size(Bond_index,1)
    Bond_index(i,4)=properties(Bond_index(i,1)).ionicradii+properties(Bond_index(i,2)).ionicradii;
    [mean_bv,std_bv,bv,bvalue]=bond_valence_data(properties(Bond_index(i,1)).type,properties(Bond_index(i,2)).type,Bond_index(i,3),Ion_1,Ion_2,R0,b,Valence_1,Valence_2,properties(Bond_index(i,1)).oxstate,properties(Bond_index(i,2)).oxstate);
    Bond_index(i,5)=bv;
end

format short
try
    diff_bond=Bond_index(:,4)-Bond_index(:,3);
    ind=find(abs(diff_bond)>0.5);
    if numel(ind)>0
        disp('Possible problems with bond between:')
        for i=1:numel(ind)
            [Bond_index(ind(i),1) Bond_index(ind(i),2)]
            properties(Bond_index(ind(i),1)).type
            properties(Bond_index(ind(i),2)).type
            Bond_index(ind(i),3:end)
        end
    end
    assignin('caller','diff_bond',diff_bond);
catch
    disp('Could not calc diff_bond')
end



assignin('caller','Tot_valence',Tot_valence);
assignin('caller','Tot_valence_oxstate',Tot_valence_oxstate);
assignin('caller','GII',GII);
assignin('caller','GII_noH',GII_noH);
assignin('caller','BondSummary',BondSummary);

assignin('caller','diff_valence',diff_valence);
assignin('caller','prop_atom',element);
try
    assignin('caller','Angle_index',Angle_index);
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','dist_matrix',dist_matrix);
    assignin('caller','diff_bond_bv',[properties.rdiffvalence]');
catch
    disp('Could not assignin Angle_index, distmatrix')
    assignin('caller','Bond_index',Bond_index);
end

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
