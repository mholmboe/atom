%% COM_atom.m
% * This function calculates the COM for certain elements
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = COM_atom(atom,Box_dim)
%
function atom = COM_atom(atom,varargin)

if nargin > 1
   % But we do not really use Box_dim do we?
   Box_dim=cell2mat(varargin(1));
end

% disp('Calculating COM')

% atom = unwrap_atom(atom,Box_dim,'xyz');

% for i=1:size(atom,2)
%     atom(i).element={atom(i).type{1}(1)};
%     atom(i).Mw=0;
% end

atom = element_atom(atom);
for i=1:size(atom,2)
    atom(i).Mw=0;
end

Element_labels=sort(unique([atom.element]));
for i=1:size(Element_labels,2)
    if strncmp(Element_labels(i),{'C'},1)
        if strncmpi(Element_labels(i),{'Cl'},2)
            [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(35.453);
        elseif strncmpi(Element_labels(i),{'Ca'},2)
            [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(40.078);
        else
            [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(12.01);
        end
    elseif strncmp(Element_labels(i),{'H'},1)
        [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(1.008);
    elseif strncmp(Element_labels(i),{'N'},1)
        [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(14.007);
    elseif strncmp(Element_labels(i),{'O'},1)
        [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(15.9994);
    elseif strncmp(Element_labels(i),{'P'},1)
        [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(30.97376);
    elseif strncmpi(Element_labels(i),{'Si'},2)
        [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(28.0855);
    elseif strncmpi(Element_labels(i),{'S'},1)
        [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(32.0650);
    elseif strncmpi(Element_labels(i),{'Al'},2)
        [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(26.981);
    elseif strncmpi(Element_labels(i),{'Mg'},2)
        [atom(find(strcmp([atom.element],Element_labels(i)))).Mw]=deal(24.305);
    else
        disp('Could not find the element this atomtype correspond to!!!')
        Element_labels(i)
    end
end


Mass_molid=sum([atom.Mw]);
for j=1:size(atom,2)
    atom(j).COM_x  = atom(j).Mw*atom(j).x/Mass_molid;
    atom(j).COM_y  = atom(j).Mw*atom(j).y/Mass_molid;
    atom(j).COM_z  = atom(j).Mw*atom(j).z/Mass_molid;
end

COM(1,1)=sum([atom.COM_x]);
COM(1,2)=sum([atom.COM_y]);
COM(1,3)=sum([atom.COM_z]);

[atom.COM_x]=deal(sum([atom.COM_x]));
[atom.COM_y]=deal(sum([atom.COM_y]));
[atom.COM_z]=deal(sum([atom.COM_z]));

assignin('caller','COM',COM);
