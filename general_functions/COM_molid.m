%% COM_molid.m
% * This function calculates the COM for a certain molid
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = COM_molid(atom,1)
%
function atom = COM_molid(atom,Box_dim,varargin)

disp('Calculating COM')
% 
% atom = unwrap_atom(atom,Box_dim,'xyz');

if nargin > 2
    MolID=varargin{1};
else
    MolID=unique([atom.molid]);
end

ind = ismember([atom.molid],MolID);

atom = element_atom(atom);

atom = mass_atom(atom);

% Element_labels=sort(unique([atom.element]))
% for i=1:size(Element_labels,2)
%     if strncmp(Element_labels(i),{'C'},1)
%         if strncmpi(Element_labels(i),{'Cl'},2)
%             [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(35.453);
%         elseif strncmpi(Element_labels(i),{'Ca'},2)
%             [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(40.078);
%         else
%             [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(12.01);
%         end
%     elseif strncmp(Element_labels(i),{'H'},1)
%         [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(1.008);
%     elseif strncmp(Element_labels(i),{'N'},1)
%         [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(14.007);
%     elseif strncmp(Element_labels(i),{'O'},1)
%         [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(15.9994);
%     elseif strncmp(Element_labels(i),{'P'},1)
%         [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(30.97376);
%     elseif strncmpi(Element_labels(i),{'Si'},2)
%         [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(28.0855);
%     elseif strncmpi(Element_labels(i),{'S'},1)
%         [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(32.0650);
%     elseif strncmpi(Element_labels(i),{'Al'},2)
%         [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(26.981);
%     elseif strncmpi(Element_labels(i),{'Mg'},2)
%         [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(24.305);
%     else
%         disp('Could not find the element this atomtype correspond to!!!')
%         Element_labels(i)
%     end
% end

COM=zeros(length(MolID),3); nCOM=1;
for i=1:numel(MolID)
    ind=find([atom.molid]==MolID(i));
    Mass_molid=sum([atom(ind).mass]);
    for j=ind
        atom(j).COM_x  = atom(j).mass*atom(j).x/Mass_molid;
        atom(j).COM_y  = atom(j).mass*atom(j).y/Mass_molid;
        atom(j).COM_z  = atom(j).mass*atom(j).z/Mass_molid;
    end
    
    COM(nCOM,1)=sum([atom(ind).COM_x]);
    COM(nCOM,2)=sum([atom(ind).COM_y]);
    COM(nCOM,3)=sum([atom(ind).COM_z]);
    
    [atom(ind).COM_x]=deal(sum([atom(ind).COM_x]));
    [atom(ind).COM_y]=deal(sum([atom(ind).COM_y]));
    [atom(ind).COM_z]=deal(sum([atom(ind).COM_z]));
    
    nCOM=nCOM+1;
end

assignin('caller','COM',COM);
