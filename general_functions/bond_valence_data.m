%% bond_valence_data.m
% * This function fetches the data and matches it to the passed atom types
% used to calculate the bond valence values according to
% * http://www.iucr.org/resources/data/datasets/bond-valence-parameters
% compiled by I. David Brown, McMaster University, Ontario, Canada
% * Data set bvparm2016.cif: 2016 version, (posted 2016-11-03)
% * atom is the atom struct
% * Box_dim is the box dimension vector
%
%% Version
% 2.081
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # [mean_bv,std_bv,bv,bvalue,ideal_dist]=bond_valence_data(ion1,ion2,R)
% # [mean_bv,std_bv,bv,bvalue,ideal_dist]=bond_valence_data(ion1,ion2,R,varargin)

function [mean_bv,std_bv,bv,bvalue]=bond_valence_data(ion1,ion2,R,varargin)

if nargin==3
    load('bond_valence_values.mat');
    valence_ion1=-100; % Dummy value
    valence_ion2=-100; % Dummy value
else
    Ion_1=varargin{1};
    Ion_2=varargin{2};
    R0=varargin{3};
    b=varargin{4};
    Valence_1=varargin{5};
    Valence_2=varargin{6};
end

if nargin>9
    valence_ion1=varargin{7};
else
    valence_ion1=-100; % Dummy value
end
if nargin>10
    valence_ion2=varargin{8};
else
    valence_ion2=-100; % Dummy value
end

if strncmpi(ion1,'Hw',2)
    ion1='H';
end

if strncmpi(ion1,'Ow',2)
    ion1='O';
end

if strncmpi(ion1,'Oh',2)
    ion1='O';
end

if strncmpi(ion2,'Hw',2)
    ion2='H';
end

if strncmpi(ion2,'Ow',2)
    ion2='O';
end

if strncmpi(ion2,'Oh',2)
    ion2='O';
end

ind1=find(strcmpi(Ion_1,ion1));
ind2=find(strcmpi(Ion_2,ion2));

if valence_ion1>-50
    valence1_ind=find(Valence_1==valence_ion1);
    ind1=intersect(ind1,valence1_ind);
end

if valence_ion2>-50
    valence2_ind=find(Valence_2==valence_ion2);
    ind2=intersect(ind2,valence2_ind);
end

ind=find(ismember(ind1,ind2));
ind=ind1(ind);

% ind
% Ion_1(ind)
% Valence_1(ind)
%
% Ion_2(ind)
% Valence_2(ind)

if numel(ind)==0
    ind1=find(ismember(Ion_2,ion1));
    ind2=find(ismember(Ion_1,ion2));
    ind=find(ismember(ind1,ind2));
    ind=ind1(ind);
end

% Ignore these H-O entries
ind(ind==636)=[];
ind(ind==637)=[];
ind(ind==638)=[];
ind(ind==639)=[];
ind(ind==640)=[];

if numel(ind)==0
    disp('Could not find any matching pair...')
    ion1
    ion2
    R
    bv=0;
    bvalue=0.37;
    bv_temp=0;
elseif (R > 1.25) && (strncmpi(ion1,'H',1) || strncmpi(ion2,'H',1))
    bv=0;
    bvalue=0.37;
    bv_temp=0;
else
    bvalue=b(ind(1));
    bv=exp((R0(ind(1))-R)/bvalue);
    for i=1:numel(ind)
        bvalue_temp=b(ind(i));
        bv_temp(i)=exp((R0(ind(i))-R)/bvalue_temp);
    end
end

mean_bv=mean(bv_temp);
std_bv=std(bv_temp);


