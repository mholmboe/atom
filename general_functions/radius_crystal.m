%% radius_crystal.m
% * For more detailed info, see the Revised Shannon radii...
%
%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # radii = radius_crystal({'O'})
% # radii = radius_crystal('O')

function radii = radius_crystal(Atom_label,varargin)

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
    try
        ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),2));
    catch
        try
            ind=find(strncmpi([Radiiproperties.Ion],Atom_label(i),1));
        catch
            ind=285;% As in O
        end
    end
    %     radii(i)=median(Radiiproperties.IonicRadii(ind
    if numel(Oxidationstate)>0
        if numel(CN)>0
            ind_cn=find(CN(i)==[Radiiproperties.CN]);
            ind_ox=find(Oxidationstate(i)==[Radiiproperties.OxState]);
            ind_cn_ox=intersect(ind_cn,ind_ox);
            ind=intersect(ind,ind_cn_ox);
            if numel(ind)>0
                radii(i,1)=Radiiproperties.CrysRadii(ind(1));
            else
               disp('Did not find this combination of ion, oxidation state and coordination number...') 
            end
        else
            ind_ox=find(Oxidationstate(i)==[Radiiproperties.OxState]);
            ind=intersect(ind,ind_ox);
            if numel(Oxidationstate)>0
                radii(i,1)=Radiiproperties.CrysRadii(ind(1));
            else
                disp('Did not find this combination of ion and oxidation state...') 
            end
        end
    else
        radii(i,1)=Radiiproperties.CrysRadii(ind(1));
    end
end



