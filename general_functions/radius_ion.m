%% radius_ion.m
% * For more detailed info, see the Revised Shannon radii...
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # radii = radius_ion({'O'})
% # radii = radius_ion('O')

function [radii,Atom_label] = radius_ion(Atom_label,varargin)

if ~iscell(Atom_label)
    Atom_label={Atom_label};
end

Oxidationstate=[];
if nargin>1
    Oxidationstate=varargin{1};
end

CN=[];
if nargin>2
    CN=varargin{2};
end

Radiiproperties=load('Revised_Shannon_radii.mat');
radii=zeros(length(Atom_label),1);
for i=1:length(Atom_label)
    
    if numel(Atom_label{i})>2
        Atom_label{i}=Atom_label{i}(1:2);
    end
    
    try
        ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),2));
    catch
        try
            ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),1));
        catch
            ind=285;% As in O
        end
    end
    ind_type=ind;
    if numel(ind_type)>1
        ind_type=ind_type(1);
    elseif numel(ind_type)==0
        try
            ind=find(strncmpi([Radiiproperties.Ion],Atom_label{i}(1),1));
            if numel(ind)==0
                ind=285;
                ind_type=285;
            else
                ind=ind(1);
                ind_type=ind;
            end
        catch
            ind=285;
            ind_type=285;
        end
    end
    %     radii(i)=median(Radiiproperties.IonicRadii(ind))';
    if numel(Oxidationstate)>0
        if numel(CN)>0
            ind_cn=find(CN(i)==[Radiiproperties.CN]);
            ind_ox=find(Oxidationstate(i)==[Radiiproperties.OxState]);
            ind_cn_ox=intersect(ind_cn,ind_ox);
            ind=intersect(ind,ind_cn_ox);
            if numel(ind)>0
                radii(i,1)=Radiiproperties.IonicRadii(ind(1));
            else
                disp('Did not find this combination of ion, oxidation state and coordination number...')
            end
        else
            ind_ox=find(Oxidationstate(i)==[Radiiproperties.OxState]);
            ind=intersect(ind,ind_ox);
            if numel(Oxidationstate)>0
                radii(i,1)=Radiiproperties.IonicRadii(ind(1));
            else
                disp('Did not find this combination of ion and oxidation state...')
            end
        end
    else
        if numel(ind)>1
            ind=ind(1);
        end
        radii(i,1)=Radiiproperties.IonicRadii(ind_type);
    end
end
Atom_label;
ind_type(1);
Atom_label=[Radiiproperties.Ion(ind_type)];
