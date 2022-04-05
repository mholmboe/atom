%% properties_atom(atom)
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
% # properties = properties_atom(atom,Box_dim)
% # properties = properties_atom(atom,Box_dim,2.5)

function properties = properties_atom(atom,Box_dim,varargin)

disp('Depreceated function... use analyze_atom(atom,Box_dim,varargin) instead')
pause(5)

% FIX multiple OxStates and other from the Shannon file

load('Revised_Shannon_radii.mat');

if nargin > 2
    rmax=varargin{1};
else
    rmax=2.25;
end

element=element_atom(atom); % Optionally, set the names ow water O and water H
element=bond_angle_atom(element,Box_dim,1.25,rmax,'more');
element=bond_valence_atom(element,Box_dim,1.25,rmax);
element=mass_atom(element);
assignin('caller','element',element);
for i=1:size(element,2)
    Ion_ind=find(strcmp([element(i).type],Ion));
    if numel(Ion_ind)==0
        Ion_ind=find(strncmpi([element(i).type],Ion,1));
    end
    CN_ind=find(numel(element(i).neigh.index)==CN);
    ind=intersect(Ion_ind,CN_ind);
    if numel(ind)==0
        ind=Ion_ind(1);
    end
    properties(i).index=element(i).index;
    properties(i).type=Ion(ind);
    properties(i).neigh=element(i).neigh;
    properties(i).bond=element(i).bond;
    properties(i).angle=element(i).angle;
    properties(i).bv=element(i).bv;
    properties(i).valence=element(i).valence;
    properties(i).rdiffvalence=element(i).Rdiff;
    properties(i).atnum=element(i).atnum;
    properties(i).mass=element(i).mass;
    properties(i).oxstate=OxState(ind);
    properties(i).ip=ZoverIR(ind);
    properties(i).cn=CN(ind);
    properties(i).crysradii=CrysRadii(ind);
    properties(i).ionicradii=IonicRadii(ind);
    properties(i).vdwradii=radius_vdw([element(i).type]);
    properties(i).elecconf=strtrim(ElecConf(ind));
    properties(i).spinstate=strtrim(SpinState(ind));
    if mod(i,100)==1
        i-1
    end
end

for i=1:size(element,2)
    if numel(properties(i).type)>1
        Ion_ind=find(strcmp([element(i).type],Ion));
        if numel(Ion_ind)==0
            Ion_ind=find(strncmpi([element(i).type],Ion,1));
        end
        CN_ind=find(numel(element(i).neigh.index)==CN);
        ind=intersect(Ion_ind,CN_ind);
        if numel(ind)==0
            ind=Ion_ind(1);
        end
        
        current_radii=[properties(i).ionicradii];
        neigh_radii=[properties([properties(i).neigh.index]).ionicradii];
        sum_radii=repmat(current_radii',numel(neigh_radii),1)+repmat(neigh_radii,numel(current_radii),1)';
        [minvalue,preox_ind]=min(abs(mean(sum_radii-mean([properties(i).neigh.dist]-properties(1).rdiffvalence))));
        
        properties(i).type=Ion(ind(preox_ind));
        properties(i).oxstate=OxState(ind(preox_ind));
        properties(i).ip=ZoverIR(ind(preox_ind));
        properties(i).cn=CN(ind(preox_ind));
        properties(i).crysradii=CrysRadii(ind(preox_ind));
        properties(i).ionicradii=IonicRadii(ind(preox_ind));
        properties(i).bv=element(i).bv;
        properties(i).valence=element(i).valence;
        properties(i).elecconf=strtrim(ElecConf(ind(preox_ind)));
        properties(i).spinstate=strtrim(SpinState(ind(preox_ind)));
    end
    if mod(i,100)==1
        i-1
    end
end


diff_valence=(abs([properties.oxstate])-[properties.valence])';
ind=find(abs(diff_valence)>0.2);
if numel(ind)>0
    disp('Possible problems with index due to multiple oxidation states,')
    disp('or atoms that are undersaturated!')
    ind
    assignin('caller','heal_ind',ind');
end

disp('Global instability index is:')
GII=(sum((abs([properties.oxstate])-[properties.valence]).^2)/size(properties,2))^0.5
if GII>0.2
    disp('GII > 0.2 --> Structure likely not super stable...');
end

load('bond_valence_values.mat');
for i=1:size(Bond_index,1)
    Bond_index(i,4)=properties(Bond_index(i,1)).ionicradii+properties(Bond_index(i,2)).ionicradii;
    [mean_bv,std_bv,bv,bvalue]=bond_valence_data(properties(Bond_index(i,1)).type,properties(Bond_index(i,2)).type,Bond_index(i,3),Ion_1,Ion_2,R0,b,Valence_1,Valence_2);
    Bond_index(i,5)=bv;
end

diff_bond=Bond_index(:,4)-Bond_index(:,3);
ind=find(abs(diff_bond)>0.05);
if numel(ind)>0
    disp('Possible problems with bond between:')
    for i=1:numel(ind)
        ind(i)
        Bond_index(ind(i),:)
    end
end

assignin('caller','GII',GII);
assignin('caller','diff_bond',diff_bond);
assignin('caller','diff_valence',diff_valence);
try
    assignin('caller','Angle_index',Angle_index);
    assignin('caller','Bond_index',Bond_index);
    assignin('caller','dist_matrix',dist_matrix);
    assignin('caller','diff_bond_bv',[properties.rdiffvalence]');
catch
    disp('Could not assignin Bond_index, Angle_index, distmatrix')
end
