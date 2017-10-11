%% lmp_atom_style_full_func.m
% * This function creates and prints the 'Atoms' properties in the LAMMPS 
% * data file.lj file according to atom style full, without image flags

function Atom_prop = lmp_atom_style_full_func(fid,Atom_label,Charge,XYZ_labels,XYZ_data) 

if size(XYZ_labels,1) > size(XYZ_data,1)
    XYZ_labels=XYZ_labels(3:end,:);
else
    XYZ_labels(:,1)=strtrim(XYZ_labels(:,1));
end

atomID = 1:size(XYZ_data,1);
molID=zeros(1,size(XYZ_data,1));

i=1;nMolID=1;
while i < size(XYZ_data,1) + 1
    if sum(ismember(Atom_label,strtrim(XYZ_labels(i)))) > 0
        if find(strncmpi('Ow',XYZ_labels(i),2))
            Atom_label_ID(i,1)=find(strncmpi(Atom_label,'Ow',2)==1);
            Atom_label_ID(i+1,1)=find(strncmpi(Atom_label,'Hw',2)==1);
            Atom_label_ID(i+2,1)=find(strncmpi(Atom_label,'Hw',2)==1);
            molID(1,i)=nMolID;
            molID(1,i+1)=nMolID;
            molID(1,i+2)=nMolID;
            i=i+2;
        else
            Atom_label_ID(i,1)=find(ismember(Atom_label,strtrim(XYZ_labels(i)))==1);
            molID(1,i)=nMolID;
        end
    else
        Atom_label_ID(i,1)=1;
    end
    nMolID=nMolID+1;
    i=i+1;
end


for i = 1:size(XYZ_data,1)
    molID(1,i);
    Atom_label_ID(i,1);
    Charge(Atom_label_ID(i,1));
    Atoms_data(i,:) = {i, molID(i), Atom_label_ID(i), Charge(Atom_label_ID(i,1)), XYZ_data(i,1),XYZ_data(i,2), XYZ_data(i,3)}; 
    fprintf(fid, '%-i\t%-i\t%-i\t%-8.6f\t%-f\t%-f\t%-f\r\n', Atoms_data{i,:});
end

Total_charge = sum(round(cell2mat(Atoms_data(:,4))*1e6)/1e6)

Atom_prop = {atomID(1:end-1), molID(1:end-1), Atom_label_ID(:,1), Charge(1,:), XYZ_data(i,1),XYZ_data(i,2), XYZ_data(i,3)};

assignin('base','atomID',atomID);
assignin('base','molID',molID);
assignin('base','Total_charge',Total_charge);





