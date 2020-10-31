%% sort_molid.m
% * This function sorts the molecular indexes in an ascending order
% 
%% Version
% 2.08
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% * sorted_molid=sort_molid(MolID)
%
function sorted_molid=sort_molid(MolID)

Tot_MolID=[];
for i=1:size(MolID,2)
   Tot_MolID=[Tot_MolID;cell2mat(MolID(i))];
end

sorted_molid={sort(Tot_MolID)};




