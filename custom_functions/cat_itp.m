function itp = cat_itp(itp1,itp2,varargin)
%% cat_itp.m
%% This function concatenates two itp structs
%% Please report bugs to michael.holmboe@umu.se

Not fininshed yet.... and not really needed
.
.
.
.
.
.
.
.
.

%% Create vars for the sections
names1 = fieldnames(itp1);
for i=1:length(names1)
    eval([strcat(names1{i},'1') '=itp1.' names1{i} ';']);
end

%% Create vars for the sections
names2 = fieldnames(itp2);
for i=1:length(names2)
    eval([strcat(names2{i},'2') '=itp2.' names2{i} ';']);
end

if size(names2,1)> size(names1,1)
    names=names2;
else
    names=names1;
end

for i=1:size(itp1.atoms.nr,1)
    atom1(i).index=itp1.atoms.nr(i);
    atom1(i).type=itp1.atoms.type(i);
    atom1(i).molid=itp1.atoms.resnr(i);
    atom1(i).resname=itp1.atoms.residue(i);
    atom1(i).fftype=itp1.atoms.atom(i);
    atom1(i).cgnr=itp1.atoms.cgnr(i);
    atom1(i).charge=itp1.atoms.charge(i);
    atom1(i).mass=itp1.atoms.mass(i);
    try
        atom1(i).typeB=itp1.atoms.typeB(i);
    catch
    end
    try
        atom1(i).chargeB=itp1.atoms.chargeB(i);
    catch
    end
    try
        atom1(i).massB=itp1.atoms.massB(i);
    catch
    end
end

for i=1:size(itp2.atoms.nr,1)
    atom2(i).index=itp2.atoms.nr(i);
    atom2(i).type=itp2.atoms.type(i);
    atom2(i).molid=itp2.atoms.resnr(i);
    atom2(i).resname=itp2.atoms.residue(i);
    atom2(i).fftype=itp2.atoms.atom(i);
    atom2(i).cgnr=itp2.atoms.cgnr(i);
    atom2(i).charge=itp2.atoms.charge(i);
    atom2(i).mass=itp2.atoms.mass(i);
    try
        atom2(i).typeB=itp2.atoms.typeB(i);
    catch
    end
    try
        atom2(i).chargeB=itp2.atoms.chargeB(i);
    catch
    end
    try
        atom2(i).massB=itp2.atoms.massB(i);
    catch
    end
end

atom=[atom1 atom2];

if exist('atoms1','var') && exist('atoms2','var')
    
end

i=1;
while i < size(names1,1)+1 | i < size(names2,1)+1
    if strcmp(names1(i))
        
    end
    i=i+1;
end


numel(fieldnames(itp))

for i=1:size(names,1)
    field=char(names(i));
    itp.(field)=Data(section_rows(i)+1:end);
end


if iscell(itp)
    % In case first struct is empty
    if size(itp,2)>1 && size(itp{1},2)==1
        for i=2:size(itp,2)
            newitp{i-1}=itp{i};
        end
        itp=newitp;
    end
    
    size_ind=zeros(size(itp,2),1);
    for i=1:size(itp,2)
        if size(itp{i},2)>0
            size_ind(i)=1;
        end
    end
    
    keepfieldnames=fieldnames(itp{1}); % Orig line
    if size(itp,2) > 1
        for i=1:size(itp,2)
            keepfieldnames=intersect(keepfieldnames,fieldnames(itp{i}));
        end
        for i=1:size(itp,2)
            rmfieldnames=setdiff(fieldnames(itp{i}),keepfieldnames);
            for j=1:numel(rmfieldnames)
                itp{i}=rmfield(itp{i},rmfieldnames(j))
            end
        end
    end
    
    itp=itp(logical(size_ind));
    itp_Tot=itp{1}; % Orig line
    if size(itp,2) > 1
        for i=2:size(itp,2)
            itp_temp=itp{i};
            if size(itp_temp,2)>0
                if numel(unique([itp_temp.molid]))==1
                    molid=num2cell(([itp_temp.molid]+[itp_Tot(end).molid]));
                    [itp_temp.molid]=deal(molid{:});
                elseif numel(unique([itp_temp.molid]))==size(itp_temp,2)
                    molid=num2cell(([1:size(itp_temp,2)]+[itp_Tot(end).molid]));
                    [itp_temp.molid]=deal(molid{:});
                end
            end
            itp_Tot=[itp_Tot itp_temp];
        end
    end
    itp=itp_Tot;
end