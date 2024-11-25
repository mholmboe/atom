function itp = modify_itp(itp,index,varargin)
%% modify_itp.m
%% This function modifies a itp struct
%% Currently it can only add values to all indexed values or remove single ones
%% Please report bugs to michael.holmboe@umu.se

%% Create vars for the sections
names = fieldnames(itp);
for i=1:length(names)
    eval([names{i} '=itp.' names{i} ';']);
end

%% If you want to add any +1 value to the atom index order
if strcmp(varargin{1},'add')
    if exist('atoms','var')
        itp.atoms.nr=itp.atoms.nr+index;
        itp.atoms.cgnr=itp.atoms.cgnr+index;
    end
    if exist('bonds','var')
        itp.bonds.ai=itp.bonds.ai+index;
        itp.bonds.aj=itp.bonds.aj+index;
    end
    if exist('angles','var')
        itp.angles.ai=itp.angles.ai+index;
        itp.angles.aj=itp.angles.aj+index;
        itp.angles.ak=itp.angles.ak+index;
    end
    if exist('pairs','var')
        itp.pairs.ai=itp.pairs.ai+index;
        itp.pairs.aj=itp.pairs.aj+index;
    end
    if exist('exclusions','var')
        itp.exclusions.ai=itp.exclusions.ai+index;
        itp.exclusions.aj=itp.exclusions.aj+index;
        itp.exclusions.ak=itp.exclusions.ak+index;
    end
    if exist('dihedrals','var')
        itp.dihedrals.ai=itp.dihedrals.ai+index;
        itp.dihedrals.aj=itp.dihedrals.aj+index;
        itp.dihedrals.ak=itp.dihedrals.ak+index;
        itp.dihedrals.al=itp.dihedrals.al+index;
    end
    if exist('impropers','var')
        itp.impropers.ai=itp.impropers.ai+index;
        itp.impropers.aj=itp.impropers.aj+index;
        itp.impropers.ak=itp.impropers.ak+index;
        itp.impropers.al=itp.impropers.al+index;
    end
end

%% If you manually deleted a single index 'index'
if strcmp(varargin{1},'delete')
    if exist('atoms','var')
        itp.atoms.nr(itp.atoms.nr>index)=itp.atoms.nr(itp.atoms.nr>index)-1;
        itp.atoms.cgnr(itp.atoms.nr>index)=itp.atoms.cgnr(itp.atoms.nr>index)-1;
    end
    if exist('bonds','var')
        itp.bonds.ai(itp.bonds.ai>index)=itp.bonds.ai(itp.bonds.ai>index)-1;
        itp.bonds.aj(itp.bonds.aj>index)=itp.bonds.aj(itp.bonds.aj>index)-1;
    end
    if exist('angles','var')
        itp.angles.ai(itp.angles.ai>index)=itp.angles.ai(itp.angles.ai>index)-1;
        itp.angles.aj(itp.angles.aj>index)=itp.angles.aj(itp.angles.aj>index)-1;
        itp.angles.ak(itp.angles.ak>index)=itp.angles.ak(itp.angles.ak>index)-1;
    end
    if exist('pairs','var')
        itp.pairs.ai(itp.pairs.ai>index)=itp.pairs.ai(itp.pairs.ai>index)-1;
        itp.pairs.aj(itp.pairs.aj>index)=itp.pairs.aj(itp.pairs.aj>index)-1;
    end
    if exist('exclusions','var')
        itp.exclusions.ai(itp.exclusions.ai>index)=itp.exclusions.ai(itp.exclusions.ai>index)-1;
        itp.exclusions.aj(itp.exclusions.aj>index)=itp.exclusions.aj(itp.exclusions.aj>index)-1;
        itp.exclusions.ak(itp.exclusions.ak>index)=itp.exclusions.ak(itp.exclusions.ak>index)-1;
    end
    if exist('dihedrals','var')
        itp.dihedrals.ai(itp.dihedrals.ai>index)=itp.dihedrals.ai(itp.dihedrals.ai>index)-1;
        itp.dihedrals.aj(itp.dihedrals.aj>index)=itp.dihedrals.aj(itp.dihedrals.aj>index)-1;
        itp.dihedrals.ak(itp.dihedrals.ak>index)=itp.dihedrals.ak(itp.dihedrals.ak>index)-1;
        itp.dihedrals.al(itp.dihedrals.al>index)=itp.dihedrals.al(itp.dihedrals.al>index)-1;
    end
    if exist('impropers','var')
        itp.impropers.ai(itp.impropers.ai>index)=itp.impropers.ai(itp.impropers.ai>index)-1;
        itp.impropers.aj(itp.impropers.aj>index)=itp.impropers.aj(itp.impropers.aj>index)-1;
        itp.impropers.ak(itp.impropers.ak>index)=itp.impropers.ak(itp.impropers.ak>index)-1;
        itp.impropers.al(itp.impropers.al>index)=itp.impropers.al(itp.impropers.al>index)-1;
    end
end

end
