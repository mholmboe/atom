%% radius_crystal.m
% * For more detailed info, see the Revised Shannon radii...
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # radii = radius_crystal({'O'})
% # radii = radius_crystal('O')

function radii = radius_crystal(Atom_label)

if ~iscell(Atom_label)
    Atom_label={Atom_label};
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
%     radii(i)=median(Radiiproperties.IonicRadii(ind))';
    radii(i,1)=Radiiproperties.CrysRadii(ind(1));
end


