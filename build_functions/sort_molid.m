%% sort_molid.m
% * This function sorts the molecular indexes in an ascending order
%
%% Version
% 3.0.1
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * sorted_molid=sort_molid(MolID)
%
%% See also
% # <replace_molid.html replace_molid>
% # <translate_molid.html translate_molid>
%
function sorted_molid=sort_molid(MolID)

Tot_MolID=[];
for i=1:size(MolID,2)
    Tot_MolID=[Tot_MolID;cell2mat(MolID(i))];
end

sorted_molid={sort(Tot_MolID)};

end




