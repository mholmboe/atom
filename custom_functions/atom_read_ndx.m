%% atom_read_ndx.m
% * This special function can read Gromacs index (.ndx) files into arrays
%
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # ndx = atom_read_ndx(filename) % Basic input arguments
% # ndx = atom_read_ndx(filename,atom) % Will generate new atom structs from the defined groups
% # ndx = atom_read_ndx(filename,atom,Box_dim,1) % Will print separate .pdb files of all groups

function ndx = atom_read_ndx(filename,varargin)

inputfile = fopen(filename, 'r');
C = textscan(inputfile, '%s', 'Delimiter', '\n');
fclose(inputfile);
D=C;

section_rows=strfind(D{1},'[');
section_rows = find(~cellfun('isempty',section_rows));
sections=D{1,1}(section_rows);

for i = 1:size(sections,1)
    sections(i)=strtrim(sections(i));
    sections(i)=strrep(sections(i),'[ ','');
    sections(i)=strrep(sections(i),' ]','');
    sections(i)=strrep(sections(i),'-','_');
end

Data=D{1,1};

i=1;
while i<size(sections,1)+1
    i;
    if i > 1
        if ismember(sections(i),sections(i+1:end))
            sections(i)
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

%% If we have duplicate sections
for i=1:size(sections,1)
    field=char(strcat(sections(i),'_ndx'));
    temp=Data(section_rows(i)+1:end);
    if i==size(sections,1)
        temp=Data(section_rows(i)+1:end);
    else
        temp=Data(section_rows(i)+1:section_rows(i+1)-1);
    end
    temp_CellArray = cellfun(@str2num,temp,'UniformOutput',false);
    NumArray=[];
    for i=1:size(temp_CellArray,1)
        NumArray=[NumArray temp_CellArray{i}];
    end
    NumArray=sort(unique(NumArray));
    ndx.(field)=NumArray;
end

%% Create vars for the sections
names = fieldnames(ndx);
for i=1:length(names)
    eval([names{i} '=ndx.' names{i} ';']);
    assignin('caller',names{i},ndx.(names{i}));
    if nargin>1
        atom=varargin{1};
        temp_name=strrep(names{i},'_ndx','');
        eval([names{i} '=ndx.' names{i} ';']);
        assignin('caller',temp_name,atom([ndx.(names{i})]));
    end
end

if nargin>3
    Box_dim=varargin{2};
    for i=1:length(names)
        temp_name=strrep(names{i},'_ndx','');
        write_atom_pdb(atom([ndx.(names{i})]),Box_dim,strcat(temp_name,'.gro'));
    end
end

