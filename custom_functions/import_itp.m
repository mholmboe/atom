function itp = import_itp(filename)
%% import_itp.m
%% This function imports a itp file Each existing section should be
%% followed by a line starting with a ';' and labels indicating the
%% parameters, something like this below...
%%
%% Written by MHolmboe
%% Please report bugs, issues
%% michael.holmboe@umu.se

% [ atomtypes ]
% ; type name atnum charge ptype v w ; v usually sigma and epsilon
%
% [ atoms ]
% ;   nr      type  resnr residue  atom   cgnr    charge      mass
%
% [ bonds ]
% ai  aj   funct  c0  c1  c2  c3
%
% [ angles ]
% ;  ai    aj    ak funct    c0    c1    c2    c3
%
% [ pairs ]
% ;  ai    aj    ak funct    c0    c1    c2    c3
%
% [ exclusions ]
% ;  ai   aj   ak   funct
%
% [ dihedrals ]
% ;  ai    aj    ak    al funct    c0    c1    c2    c3    c4    c5
%
% [ impropers ]
% ;  ai   aj   ak   al  funct   c0   c1   c2   c3   c4   c5

inputfile = fopen(filename, 'r');
C = textscan(inputfile, '%s', 'Delimiter', '\n');
fclose(inputfile);
D=C;
nondatarows1=find(strncmp(C{1},';',1));
nondatarows2=find(strncmp(C{1},{''},1));
nondatarows3=find(strncmp(C{1},{'#'},1));
D{1,1}([nondatarows1;nondatarows2;nondatarows3])=[];

section_rows=strfind(D{1},'[');
section_rows = find(~cellfun('isempty',section_rows));
sections=D{1,1}(section_rows);

Data=D{1,1};

i=1;
while i<size(sections,1)+1
    i;
    sections(i)=strtrim(sections(i));
    if i > 1
        if ismember(sections(i),sections(1:i-1))
            sections(i);
            sections(i)=[];
            Data(section_rows(i))=[];
            section_rows(i)=[];
            i=i-1;
            try
                section_rows(i+1:end)=section_rows(i+1:end)-1;
            catch
            end
        end
    end
    i=i+1;
end

for i = 1:size(sections,1)
    sections(i)=strtrim(sections(i));
    sections(i)=strrep(sections(i),'[ ','');
    sections(i)=strrep(sections(i),' ]','');
end

%% If we have duplicate sections
for i=1:size(sections,1)
    field=char(sections(i))
    field(regexp(field,';'):end)=[];
    field(regexp(field,' '):end)=[];
    if i==size(sections,1)
        itp.(field)=Data(section_rows(i)+1:end);
    else
        itp.(field)=Data(section_rows(i)+1:section_rows(i+1)-1);
    end
end

%% Create vars for the sections
names = fieldnames(itp);
for i=1:length(names)
    eval([names{i} '=itp.' names{i} ';']);
end

try
    %% Parse the [ moleculetype ] section
    if exist('moleculetype','var')
        clearvars moleculetype
        atoms_cell = regexp(itp.moleculetype,'\s+', 'split');
        atoms_cell = vertcat(atoms_cell{:});
        atoms_section_labels = {'moleculetype' 'nrexcl'};
        for i=1:size(atoms_cell,2)
            field=char(atoms_section_labels(i));
            %       if ismember(atoms_section_labels(i),{'moleculetype' 'nrexcl'})
            %         bonds.(field)=cellfun(@str2double,atoms_cell(:,i));
            %         else
            moleculetype.(field)=atoms_cell(:,i);
            %         end
        end
        itp.moleculetype=moleculetype;
    end
catch
    disp('Found no moleculetype section')
end

%% Parse the [ atoms ] section
if exist('atoms','var')
    clearvars atoms
    atoms_cell = regexp(itp.atoms,'\s+', 'split');
    atoms_cell = vertcat(atoms_cell{:});
    atoms_section_labels = {'nr' 'type' 'resnr' 'residue' 'atom' 'cgnr' 'charge' 'mass' 'typeB' 'chargeB' 'massB'  'comment'};
    for i=1:min([12 size(atoms_cell,2)])
        field=char(atoms_section_labels(i));
        if ismember(atoms_section_labels(i),{'nr' 'resnr' 'cgnr' 'charge' 'mass'})
            atoms.(field)=cellfun(@str2double,atoms_cell(:,i));
        else
            atoms.(field)=atoms_cell(:,i);
        end
    end
    itp.atoms=atoms;
end

try
    %% Parse the [ atomtypes ] section
    if exist('atomtypes','var')
        try
            clearvars atomtypes
            atoms_cell = regexp(itp.atomtypes,'\s+', 'split');
            atoms_cell = vertcat(atoms_cell{:});
            if numel(atoms_cell)>0
                if sum(strcmp([atoms_cell(1,:)],';'))>0
                    atoms_cell = atoms_cell(:,1:find(strcmp([atoms_cell(1,:)],';'))-1);
                end
                atoms_cell = sortrows(atoms_cell,1);
                atoms_section_labels = {'type' 'name' 'atnum' 'charge' 'ptype' 'v' 'w'};
                for i=1:min([7 size(atoms_cell,2)])
                    field=char(atoms_section_labels(i));
                    if ismember(atoms_section_labels(i),{'atnum' 'charge' 'v' 'w'})
                        atomtypes.(field)=cellfun(@str2double,atoms_cell(:,i));
                    else
                        atomtypes.(field)=atoms_cell(:,i);
                    end

                end
            end
            itp.atomtypes=atomtypes;
            for i=1:size(itp.atoms.type,1)
                itp.atoms.typeind(i)=find(strcmp([itp.atoms.type(i)],[itp.atomtypes.type]));
            end
        catch
            clearvars atomtypes
            atoms_cell = regexp(itp.atomtypes,'\s+', 'split');

            colmin=100;
            for i=1:size(atoms_cell,1)
                colmin=min([size(atoms_cell{i},2) colmin]);
            end
            for i=1:size(atoms_cell,1)
                atoms_cell{i}=atoms_cell{i}(1:colmin);
            end
            atoms_cell = vertcat(atoms_cell{:});
            if numel(atoms_cell)>0
                if sum(strcmp([atoms_cell(1,:)],';'))>0
                    atoms_cell = atoms_cell(:,1:find(strcmp([atoms_cell(1,:)],';'))-1);
                end
                atoms_cell = sortrows(atoms_cell,1);
                atoms_section_labels = {'name' 'atnum' 'mass' 'charge' 'ptype' 'sigma' 'epsilon'};
                for i=1:min([7 size(atoms_cell,2)])
                    field=char(atoms_section_labels(i));
                    if ismember(atoms_section_labels(i),{'atnum' 'charge' 'sigma' 'epsilon'})
                        atomtypes.(field)=cellfun(@str2double,atoms_cell(:,i));
                    else
                        atomtypes.(field)=atoms_cell(:,i);
                    end

                end
            end
            itp.atomtypes=atomtypes;
            itp.atomtypes.type=itp.atomtypes.name;
            for i=1:size(itp.atoms.type,1)
                itp.atoms.typeind(i)=find(strcmp([itp.atoms.type(i)],[itp.atomtypes.type]));
            end
        end
    end
catch
    disp('Found no atomtypes section')
end

%% Parse the [ bonds ] section
if exist('bonds','var')
    clearvars bonds
    atoms_cell = regexp(itp.bonds,'\s+', 'split');
    atoms_cell = vertcat(atoms_cell{:});
    if numel(atoms_cell)>0
        if sum(strcmp([atoms_cell(1,:)],';'))>0
            atoms_cell = atoms_cell(:,1:find(strcmp([atoms_cell(1,:)],';'))-1);
        end
        atoms_section_labels = {'ai' 'aj' 'funct' 'c0' 'c1' 'c2' 'c3'};
        for i=1:min([5 size(atoms_cell,2)])
            field=char(atoms_section_labels(i));
            if ismember(atoms_section_labels(i),{'ai' 'aj' 'funct' 'c0' 'c1' 'c2' 'c3'})
                bonds.(field)=cellfun(@str2double,atoms_cell(:,i));
                [bonds.(field)(isnan(bonds.(field)))]=deal(0);
            else
                bonds.(field)=atoms_cell(:,i:end);
                [bonds.(field)(isnan(bonds.(field)))]=deal(0);
            end

        end
    end
    if exist('bonds','var')
        itp.bonds=bonds;
    end
end

%% Parse the [ angles ] section
if exist('angles','var')
    clearvars angles
    atoms_cell = regexp(itp.angles,'\s+', 'split');

    %% To avoid vertcat problems
    mincol=1000;
    for s=1:size(atoms_cell,1)
        temp=size(atoms_cell{s},2);
        if temp < mincol
            mincol=temp;
        end
    end
    for s=1:size(atoms_cell,1)
        atoms_cell{s}(1:mincol);
    end

    atoms_cell = vertcat(atoms_cell{:});
    if numel(atoms_cell)>0
        if sum(strcmp([atoms_cell(1,:)],';'))>0
            atoms_cell = atoms_cell(:,1:find(strcmp([atoms_cell(1,:)],';'))-1);
        end

        atoms_section_labels = {'ai' 'aj' 'ak' 'funct' 'c0' 'c1' 'c2' 'c3'};
        for i=1:min([6 size(atoms_cell,2)])
            field=char(atoms_section_labels(i));
            if ismember(atoms_section_labels(i),{'ai' 'aj' 'ak' 'funct' 'c0' 'c1' 'c2' 'c3'})
                angles.(field)=cellfun(@str2double,atoms_cell(:,i));
                [angles.(field)(isnan(angles.(field)))]=deal(0);
            else
                angles.(field)=atoms_cell(:,i);
                [angles.(field)(isnan(angles.(field)))]=deal(0);
            end
        end
    end
    if exist('angles','var')
        itp.angles=angles;
    end
end

try
    %% Parse the [ pairs ] section
    if exist('pairs','var')
        clearvars pairs
        atoms_cell = regexp(itp.pairs,'\s+', 'split');
        atoms_cell = vertcat(atoms_cell{:});
        if numel(atoms_cell)>0
            if sum(strcmp([atoms_cell(1,:)],';'))>0
                atoms_cell = atoms_cell(:,1:find(strcmp([atoms_cell(1,:)],';'))-1);
            end

            atoms_section_labels = {'ai' 'aj' 'funct' 'c0' 'c1' 'c2' 'c3'};
            for i=1:min([7 size(atoms_cell,2)])
                field=char(atoms_section_labels(i));
                if ismember(atoms_section_labels(i),{'ai' 'aj' 'funct' 'c0' 'c1' 'c2' 'c3'})
                    pairs.(field)=cellfun(@str2double,atoms_cell(:,i));
                    %             [pairs.(field)(isnan(pairs.(field)))]=deal(0);
                else
                    pairs.(field)=atoms_cell(:,i);
                    %             [pairs.(field)(isnan(pairs.(field)))]=deal(0);
                end
            end
        end
        if exist('pairs','var')
            itp.pairs=pairs;
        end

    end
catch
end

try
    %% Parse the [ exclusions ] section
    if exist('exclusions','var')
        clearvars exclusions
        atoms_cell = regexp(itp.exclusions,'\s+', 'split');
        atoms_cell = vertcat(atoms_cell{:});
        if numel(atoms_cell)>0
            if sum(strcmp([atoms_cell(1,:)],';'))>0
                atoms_cell = atoms_cell(:,1:find(strcmp([atoms_cell(1,:)],';'))-1);
            end

            atoms_section_labels = {'ai' 'aj' 'ak' 'funct'};
            for i=1:min([4 size(atoms_cell,2)])
                field=char(atoms_section_labels(i));
                if ismember(atoms_section_labels(i),{'ai' 'aj' 'ak' 'funct'})
                    exclusions.(field)=cellfun(@str2double,atoms_cell(:,i));
                else
                    exclusions.(field)=atoms_cell(:,i);
                end
            end
        end
        if exist('exclusions','var')
            itp.exclusions=exclusions;
        end
    end
catch
end

try
    %% Parse the [ dihedrals ] section
    if exist('dihedrals','var')
        clearvars dihedrals
        atoms_cell = regexp(itp.dihedrals,'\s+', 'split');
        atoms_cell = vertcat(atoms_cell{:});
        if numel(atoms_cell)>0
            if sum(strcmp([atoms_cell(1,:)],';'))>0
                atoms_cell = atoms_cell(:,1:find(strcmp([atoms_cell(1,:)],';'))-1);
            end

            atoms_section_labels = {'ai' 'aj' 'ak'  'al'  'funct' 'c0' 'c1' 'c2' 'c3'  'c4' 'c5'};
            for i=1:min([11 size(atoms_cell,2)])
                field=char(atoms_section_labels(i));
                if ismember(atoms_section_labels(i),{'ai' 'aj' 'ak'  'al'  'funct' 'c0' 'c1' 'c2' 'c3'  'c4' 'c5'})
                    dihedrals.(field)=cellfun(@str2double,atoms_cell(:,i));
                else
                    dihedrals.(field)=char(atoms_cell(:,i));
                end
            end
        end
        if exist('dihedrals','var')
            itp.dihedrals=dihedrals;
        end
    end
catch
end

try
    %% Parse the [ impropers ] section
    if exist('impropers','var')
        clearvars impropers
        atoms_cell = regexp(itp.impropers,'\s+', 'split');
        atoms_cell = vertcat(atoms_cell{:});
        if numel(atoms_cell)>0
            if sum(strcmp([atoms_cell(1,:)],';'))>0
                atoms_cell = atoms_cell(:,1:find(strcmp([atoms_cell(1,:)],';'))-1);
            end

            atoms_section_labels = {'ai' 'aj' 'ak'  'al'  'funct' 'c0' 'c1' 'c2' 'c3'  'c4' 'c5'};
            for i=1:min([11 size(atoms_cell,2)])
                field=char(atoms_section_labels(i));
                if ismember(atoms_section_labels(i),{'ai' 'aj' 'ak'  'al'  'funct' 'c0' 'c1' 'c2' 'c3'  'c4' 'c5'})
                    impropers.(field)=cellfun(@str2double,atoms_cell(:,i));
                else
                    impropers.(field)=atoms_cell(:,i);
                end
            end
        end
        if exist('impropers','var')
            itp.impropers=impropers;
        end
    end
catch
end

try
    %% Parse the [ position_restrains ] section, starting three lines above the first [ position_restraints ] instance
    if exist('position_restraints','var')

        inputfile = fopen(filename, 'r');
        PR = textscan(inputfile, '%s', 'Delimiter', '\n');
        fclose(inputfile);

        section_rows=strfind(PR{1},'position_restraints');
        section_rows = find(~cellfun('isempty',section_rows));
        sections=PR{1,1}(section_rows);
        End_Data=PR{1,1};
        End_Data=char(End_Data(section_rows(1)-3:size(PR{1,1},1),1));
        itp.enddata=End_Data;

    end
catch
end

% assignin('caller','itp',itp);

end


