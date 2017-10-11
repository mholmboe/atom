%% write_atom_mol2.m
% * This function writes an mol2 file from the atom struct
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * write_atom_mol2(atom,Box_dim,filename_out,1.25,1.25,'clayff','spce')

function write_atom_mol2(atom,Box_dim,filename_out,varargin)

filename_out=strcat(filename_out,'.mol2');

if nargin > 3
    short_r=cell2mat(varargin(1));
    long_r=cell2mat(varargin(2));
else
    short_r=1.25;
    long_r=2.25;
end

if nargin>5;
    ffname=varargin(3);
    if nargin>6;
        watermodel=varargin(4);
    else
        disp('Unknown watermodel, will try SPC/E')
        watermodel='SPC/E';
    end
    if strncmpi(ffname,'clayff',5);
        clayff_param(sort(unique([atom.type])),watermodel);
        Total_charge = check_clayff_charge(atom)
    elseif strcmpi(ffname,'interface');
        interface_param(sort(unique([atom.type])),watermodel);
        Total_charge = check_interface_charge(atom)
    elseif strcmpi(ffname,'interface15');
        interface15_param(sort(unique([atom.type])),watermodel);
        Total_charge = check_interface_charge(atom)
    else
        disp('Unknown forcefield, will try clayff')
        clayff_param(sort(unique([atom.type])),watermodel);
        Total_charge = check_clayff_charge(atom)
    end
end

% This file could write out the Box_dim as a comment!
atom=bond_angle_atom(atom,Box_dim,short_r,long_r);
residuename='MOL';
%filename_out='MMT_matlab.mol2';
nAtoms=size(atom,2);
nBonds=size(Bond_index,1);

Atom_label=sort(unique([atom.type]'));
Element=[atom.type]'; Atom_label_ID=zeros(size(Element,1),1);
% clayff_param(unique([atom(:).type]),'SPC/E');
for i=1:size(Element,1)
    Atom_label_ID(i,1)=find(ismember(Atom_label,strtrim(atom(i).type))==1);
    atom(i).charge=Charge(Atom_label_ID(i,1));
    Element(i)={Element{i}(1)};
    if strcmp(Element(i),{'A'})
        Element(i)={'Al'};
    elseif strcmp(Element(i),{'M'})
        Element(i)={'Mg'};
    elseif strcmp(Element(i),{'S'})
        Element(i)={'Si'};
    elseif strcmp(Element(i),{'F'})
        Element(i)={'Fe'};
    elseif strcmp(Element(i),{'H'})
        Element(i)={'H'};
    elseif strcmp(Element(i),{'O'})
        Element(i)={'O'};
    elseif strcmp(Element(i),{'N'})
        Element(i)={'Na'};
    elseif strcmp(Element(i),{'C'})
        Element(i)={'Ca'};
    elseif strcmp(Element(i),{'K'})
        Element(i)={'K'};
    else
        Element(i)=Element(i);
    end
end

nAtoms=size([atom.x],2)
Atom_section=cell(nAtoms,10);
fid = fopen(filename_out, 'wt');
fprintf(fid, '@<TRIPOS>MOLECULE\n');
fprintf(fid, 'MMT\n');
fprintf(fid, '%5i%6i%6i%6i%6i\n',[nAtoms nBonds 1  0 1]);
fprintf(fid, 'SMALL\n');
fprintf(fid, 'USER_CHARGES\n');
fprintf(fid, '@<TRIPOS>ATOM\n');

for i = 1:nAtoms
    Atom_section(1:10) = [atom(i).index atom(i).type atom(i).x atom(i).y atom(i).z Element(i) 1 residuename atom(i).charge '****'];
    fprintf(fid, '%3i%3s%18.6f%12.6f%12.6f%3s%10i%4s%14.4f%5s\n', Atom_section{1:10});
end
fprintf(fid, '@<TRIPOS>BOND\n');
for i = 1:size(Bond_index,1)
    Bonds_section = [i Bond_index(i,1) Bond_index(i,2) 1];
    fprintf(fid, '%5i%6i%6i%2i\n', Bonds_section);
end
fprintf(fid, '@<TRIPOS>SUBSTRUCTURE\n');
fprintf(fid, '      1  RES              1 ****               0 ****  ****');
fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

fclose(fid);
disp('.mol2 structure file written')
