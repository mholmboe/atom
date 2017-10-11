%% import_xyz.m - This function imports an .xyz file. Atom types should be made of letters, not numbers...
function import_xyz(XYZ_filename)
%% Try the import_atom_xyz function instead...


%[A,delimiterOut]=importdata(XYZ_filename);
newData1 = importdata(XYZ_filename,' ', 2); % \t

if iscell(newData1(1)) == 0 && size(newData1.data,2) == 3;
    disp('Letters in xyz-file, space separated')
    newData1 = importdata(XYZ_filename,' ', 2);
    vars = fieldnames(newData1);
    for i = 1:length(vars)
        assignin('base', vars{i}, newData1.(vars{i}));
    end
else
    disp('Letters in xyz-file, tabb separated')
    newData1 = importdata(XYZ_filename, '\t', 2);
    vars = fieldnames(newData1);
    for i = 1:length(vars)
        assignin('base', vars{i}, newData1.(vars{i}));
    end
end

XYZ_data = newData1.(vars{1});
XYZ_labels = newData1.(vars{2});

if size(XYZ_labels,1) > size(XYZ_data,1);
    XYZ_labels=XYZ_labels(3:end,:);
end

assignin('caller', 'XYZ_data', XYZ_data);
assignin('caller', 'XYZ_labels', XYZ_labels);

