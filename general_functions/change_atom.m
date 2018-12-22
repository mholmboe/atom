clear all;
format short;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename='new_8x7MMT_Skipper.gro'; Title=filename;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot = 1; % 1/0 for yes/no
writepdb  =  0; % 1/0 for yes/no
writemol2 =  0; % 1/0 for yes/no
writexyz  =  0; % 1/0 for yes/no
writegro  =  1; % 1/0 for yes/no
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if regexp(filename,'.pdb') > 1; import_atom(filename); filename_out=regexprep(filename,'.pdb',' '); disp('Found .pdb file'); end
if regexp(filename,'.mol2') > 1; filename_out=regexprep(filename,'.mol2',' ');disp('Found .mol2 file'); end
if regexp(filename,'.gro') > 1; import_atom_gro(filename); filename_out=regexprep(filename,'.gro',' '); disp('Found .gro file'); end
if regexp(filename,'.xyz') > 1; import_atom_xyz(filename); filename_out=regexprep(filename,'.xyz',' '); disp('Found .xyz file');
    Box_dim = [42.210    120  37.800000	0.000000	0.000000	0.000000]; % optional, [x_hi y_hi z_hi xy xz yz] in Ångströms
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename_out=strcat('new_',filename_out);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atom = ave_atom_func(atom);
% atom = median_atom_func(atom);
atom = wrap_molid_func(atom,Box_dim)
% atom = unwrap_atom_func(atom,Box_dim,'z');

atom = atom_update_func(atom);
% atom = translate_atom_func(atom,[0 0 0],'MMT')
% atom = wrap_atom_func(atom,Box_dim);
% atom = center_atom_func(atom,[Box_dim(1) Box_dim(2) Box_dim(3) 0 0 0],'all','z')
% atom = scale_atom_func(atom,[.64/.575 .64/.575 1],Box_dim,'all')
% atom = remove_molid_func(atom,MolID)
% atom = remove_residues_func(atom,'OW',0,120,'z')
% atom = remove_type_func(atom,typescell)
% atom = remove_resname_func(atom,{'SOL' 'NA' 'CL'});
% atom = COM_atom_func(atom,Box_dim,72);
% atom = translate_molid_func(atom,[0 0 Box_dim(3)/2-COM_molid(end,3)],72);
% atom = ave_atom_func(atom);
% atom = median_atom_func(atom);
% atom = wrap_molid_func(atom,Box_dim)
% atom = unwrap_atom_func(atom,Box_dim,'z');
% atom = atom_update_func(atom,MolID);
% atom = translate_atom_func(atom,[0 0 0],'MMT')
% atom = wrap_atom_func(atom,Box_dim);
% atom = center_atom_func(atom,[4.567 4.567 120 0 0 0],'all','z')
% atom = scale_atom_func(atom,[.64/.575 .64/.575 1],Box_dim,'all')
% atom = remove_molid_func(atom,MolID)
% atom = remove_residues_func(atom,'OW',0,120,'z')
% atom = remove_type_func(atom,typescell)
% atom = remove_resname_func(atom,{'SOL' 'NA' 'CL'});
% atom = COM_atom_func(atom,Box_dim,72);
% atom = translate_molid_func(atom,[0 0 -50],2);
% atom = duplicate_atom_func(atom,1); % atom,molid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script translates a certain residue

% ind_resname=find(strcmp([atom.resname],'TCHH'));
% ave_atom_func(atom(ind_resname));
% translate_atom_func(atom,-[[26.9567500000000,18.4087500000000,7.01150000000000]],'TCHH');
% translate_atom_func(atom,[Box_dim(1)/2,Box_dim(1)/2,10],'TCHH');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script can be used to sort all atoms based on median z of each residue
%
% atom_ind_lo_z=find([atom.med_z]<=Box_dim(3)/2);
% atom_ind_hi_z=find([atom.med_z]>Box_dim(3)/2);
% atom=atom(atom_ind_hi_z);
% atom=atom_update_func(atom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script is used to remove residue with resid/molids

% atom=remove_resname_func(atom,{'SOL' 'NA' 'CL'});

% ind_rm = setdiff(unique([atom.molid]),[1:64 219 258 293 261 198 273 122 319]);
% atom = remove_molid_func(atom,ind_rm);
% atom = median_atom_func(atom);
% atom = unwrap_atom_func(atom,Box_dim,'xyz');
% atom = translate_atom_func(atom,[0 0 Box_dim(3)*.6],'DLiPC');
% atom = translate_molid_func(atom,[0 0 (160-530)],1:32);
% atom = condense_atom_func(atom,Box_dim,3); % 5.290000 9.920000 18.139999 90.000000 90.000000 90.000000
% atom = unwrap_atom_func(atom,Box_dim,'xyz');
% atom = center_atom_func(atom,Box_dim,'DLiPC','z')
% atom = remove_molid_func(atom,[2 3 4 70 71 72]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script is used to remove residue with atomytypes...

%atom=remove_type_func(atom,{'H12S' 'H13S' 'H12Y' 'H13Y'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script is used to solvate with SPC water, optionally remove water somewhere and then
%% dilute the number of h2o to a certain number
% nRes1=64;nRes2=8;
% nH2O=30*(nRes1+nRes2);
% Ion1='NA';nIon1=ceil(0.15*nH2O/55.3);
% Ion2='CL';nIon2=ceil(0.15*nH2O/55.3);
%
% atom=center_atom_func(atom,[Box_dim(1) Box_dim(2) 90],'DLiPC','z');write_atom_gro(atom,Box_dim,'temp');
% lo_limit=Box_dim(3)/2-25;
% hi_limit=Box_dim(3)/2+25;
% str=strcat('genbox -cp',{' '},'temp',{' '},'-cs spc216.gro',{' '},'-o temp.gro &>/dev/null');
% system(strcat(char({PATH2GMX}),char(str)));
% import_atom_gro('temp.gro');
% atom=remove_residues_func(atom,'OW',lo_limit,hi_limit,'z');
% atom=dilute_water_func(atom,nH2O);
% nSOL=size(find(strcmp([atom.type],'OW')),2)
% atom=add_ions_func(atom,'OW','z',Ion1,nIon1,Ion2,nIon2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script is used to remove residues
%% in the simulation box between zlo and zhi
% %
% atomname='OW';
% xlo=10;ylo=10;zlo=10;
% xhi=50;yhi=50;zhi=50;
%
% nAtoms=size([atom.x],2);
% ind=find(strcmp([atom.type],atomname));
% ind_lo=find([atom.x]>xlo);ind_lo=intersect(ind_lo,ind);
% ind_hi=find([atom.x]<xhi);ind_hi=intersect(ind_hi,ind);
% ind=intersect(ind_lo,ind_hi);
% if strcmp(atomname,'OW')
%     ind_HW1=ind+1;ind_HW2=ind+2;
%     ind=sort([ind ind_HW1 ind_HW2]);
% end
% change_ind=sort(setdiff([1:nAtoms],ind));
% atom=atom(change_ind);
%
% nAtoms=size([atom.x],2);
% ind=find(strcmp([atom.type],atomname));
% ind_lo=find([atom.y]>ylo);ind_lo=intersect(ind_lo,ind);
% ind_hi=find([atom.y]<yhi);ind_hi=intersect(ind_hi,ind);
% ind=intersect(ind_lo,ind_hi);
% if strcmp(atomname,'OW')
%     ind_HW1=ind+1;ind_HW2=ind+2;
%     ind=sort([ind ind_HW1 ind_HW2]);
% end
% change_ind=sort(setdiff([1:nAtoms],ind));
% atom=atom(change_ind);
%
% nAtoms=size([atom.x],2);
% ind=find(strcmp([atom.type],atomname));
% ind_lo=find([atom.z]>zlo);ind_lo=intersect(ind_lo,ind);
% ind_hi=find([atom.z]<zhi);ind_hi=intersect(ind_hi,ind);
% ind=intersect(ind_lo,ind_hi);
% if strcmp(atomname,'OW')
%     ind_HW1=ind+1;ind_HW2=ind+2;
%     ind=sort([ind ind_HW1 ind_HW2]);
% end
% change_ind=sort(setdiff([1:nAtoms],ind));
% atom=atom(change_ind);
%
% center_atom_func(atom,[10 10 10 0 0 0],'all','y');
% translate_atom_func(atom,-Box_dim/2,'all')
%
% nAtoms=size(atom,2);
% size(find(strcmp([atom.resname],atomname)),2);
% nSOL=size(find(strcmp([atom.type],'OW')),2)
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% What does this script do??

% nresname=0;j=1;
% for i=1:nAtoms;
%     if i == 1
%         nresname=nresname+1;
%         ResNameTable(j,1)=atom(i).resname;
%     elseif i > 1 && strcmp(atom(i).resname,atom(i-1).resname)==0;
%         ResNameTable(j,2)={nresname};
%         j=j+1;
%         ResNameTable(j,1)=atom(i).resname;
%         nresname=1;
%     elseif i > 1
%         nresname=nresname+1;
%     else
%     end
% end
% ResNameTable(j,2)={nresname};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script removes atoms below certain z

% ind_keep=find([atom.z]>25);
% atom=atom(ind_keep);
% ind_keep=find([atom.z]<45);
% atom=atom(ind_keep);
% nAtoms=size(atom,2);
% ind_Si=find(strncmp([atom.type],'Si',2));
% z_shift=num2cell([atom.z]-mean([atom(ind_Si).z])); [atom(:).z]=deal(z_shift{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script sets the residue length with regard to nAtoms in resname 1
% nAtoms=size(atom,2);
% natomsinresidue=length(find([atom(:).molid]==1));
% for i=1:size(atom,2)
%     atom(i).molid=ceil(i/natomsinresidue);
%     %% Test this: atom(i).molid=mod(i,54)(+1?);
%     atom(i).index=mod(i,100000);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script is used to translate residues
%% in the simulation box
% for i=1:max([atom(:).index])
%     atom(i).molid
%     if atom(i).molid==129;
%         i
%     atom(i).x=atom(i).x;
%     atom(i).y=atom(i).y;
%     atom(i).z=atom(i).z-5;
%     end
% end
% i

% Box_dim=[4.28552   4.28552  11.20336]
%
% Residue= 'TCH'; % resname to be translated
% x_ave=mean([atom.x]);
% y_ave=mean([atom.y]);
% z_ave=mean([atom.z]);
% % x_shift=0;
% % y_shift=0;
% % z_shift=0;

%x_shift=num2cell([[atom.x]-mean([atom.x]-Box_dim(1)/2)]'); [atom(:).x]=deal(x_shift{:});
%x_shift=num2cell([[atom.x]-Box_dim(1)/2]'); [atom(:).x]=deal(x_shift{:});
% ind_Si=find(strncmp([atom.type],'Si',2));
%y_shift=num2cell([atom.y]-(mean([atom(ind_Si).y]-Box_dim(2)/2))); [atom(:).y]=deal(y_shift{:});
% z_shift=num2cell([atom.z]-mean([atom(ind_Si).z])); [atom(:).z]=deal(z_shift{:});
%z_shift=num2cell([[atom.z]-min([atom.z])]'); [atom(:).z]=deal(z_shift{:});
%z_shift=num2cell([[atom.z]-Box_dim(3)/6]'); [atom(:).z]=deal(z_shift{:});
%
% S=num2cell(1:39);
% [atom(:).index]=deal(S{:});

% i=1; Index=zeros(size(dataArray{1,1}(:,1),1),1);
% while i<size(Index,1)
%     if strcmp(strtrim(Resname{i}),'TCH');% == 1 & (Z_coord(i) > Zlo) & (Z_coord(i) < Zhi);
%         i
%         dataArray{:,5}(i)=dataArray{:,5}(i)+x_shift;
%         dataArray{:,6}(i)=dataArray{:,6}(i)+y_shift;
%         dataArray{:,7}(i)=dataArray{:,7}(i)+z_shift;
%     end
%     i=i+1;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script is used to remove all but 'Resname' from the stated
% Resname1='DLiPC'
% %Resname2='TCH'
% %Resname3='SOL'
% ind1=find(strcmp([atom.resname],Resname1));
% %ind2=find(strcmp([atom.resname],Resname2));
% %ind3=find(strcmp([atom.resname],Resname3));
% ind=[ind1];
% atom=atom(ind);MolID=MolID(ind);nAtoms=size(atom,2);
% nAtoms=size([atom.x],2);
%
% atom=atom_update_func(atom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I think this script is used to remove all but 'Resname' from and moliID's from the stated conf

% ind_Resname=find(strcmp([atom.resname],Resname));
% ind_molid=molid=[atom(ind_Resname).molid]; molid=unique(molid);
%
% keep_molid=[134 137 167 168 195 199 202 228];keep_ind=[];
% for i=1:length(keep_molid);
%     find([atom.molid]==keep_molid(i));
%     keep_ind=[keep_ind find([atom.molid]==keep_molid(i))];
% end
% keep_ind = sort(keep_ind);
% atom=atom(keep_ind);MolID=MolID(keep_ind);nAtoms=size(atom,2);
% atom=atom_update_func(atom,MolID);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script is used to move water
%% from somewhere in Z to somewhere else in Z

% atomname='OW';
% zlo=37;
% zhi=63;
% z_translate=30;
% zmiddle=Box_dim(3)/2;
% nAtoms=size([atom.x],2);
% ind=find(strcmp([atom.type],atomname));
% ind_lo=find([atom.z]>zlo);ind_lo=intersect(ind_lo,ind);
% ind_hi=find([atom.z]<zmiddle);ind_hi=intersect(ind_hi,ind);
% ind=intersect(ind_lo,ind_hi);
% if strcmp(atomname,'OW')
%     ind_HW1=ind+1;ind_HW2=ind+2;
%     ind=sort([ind ind_HW1 ind_HW2]);
% end
% z_shift=num2cell([[atom(ind).z]-z_translate]'); [atom(ind).z]=deal(z_shift{:});
% disp('Moving some SOL down')
% size(ind,2)/3
%
% % And repeat...
% ind=find(strcmp([atom.type],atomname));
% ind_lo=find([atom.z]>zmiddle);ind_lo=intersect(ind_lo,ind);
% ind_hi=find([atom.z]<zhi);ind_hi=intersect(ind_hi,ind);
% ind=intersect(ind_lo,ind_hi);
% if strcmp(atomname,'OW')
%     ind_HW1=ind+1;ind_HW2=ind+2;
%     ind=sort([ind ind_HW1 ind_HW2]);
% end
% z_shift=num2cell([[atom(ind).z]+z_translate]'); [atom(ind).z]=deal(z_shift{:});
% disp('Moving some SOL upp')
% size(ind,2)/3

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script adds ions on top of water
%% in the simulation box
%
% Ion1='NA'
% Ion2='CL'
% replace_atom='OW'
% nIon1=129+38;
% nIon2=38;
%
% nIon1=nIon1+nIon2;
% ind_rep=find(strcmp([atom.type],replace_atom));
% ind_add_Ion1=ind_rep(randperm(length(ind_rep),nIon1));
% Ion1_atom=atom(ind_add_Ion1);
%
% XCoords=num2cell([[Ion1_atom.x]-0.5]');
% YCoords=num2cell([[Ion1_atom.y]-0.5]');
% ZCoords=num2cell([[Ion1_atom.z]-0.5]');
%
% [Ion1_atom(:).x]=deal(XCoords{:});
% [Ion1_atom(:).y]=deal(YCoords{:});
% [Ion1_atom(:).z]=deal(ZCoords{:});
%
% [Ion1_atom(:).resname]=deal({Ion1});
% [Ion1_atom(:).type]=deal({Ion1});
% [Ion1_atom(:).fftype]=deal({Ion1});
%
% ind_Ion1=find(strcmp([Ion1_atom.type],Ion1));
% ind_add_Ion2=ind_Ion1(randperm(length(ind_Ion1),nIon2));
% Ion2_atom=Ion1_atom(ind_add_Ion2);
% Ion1_atom(ind_add_Ion2)=[];
%
% [Ion2_atom(:).resname]=deal({Ion2});
% [Ion2_atom(:).type]=deal({Ion2});
% [Ion2_atom(:).fftype]=deal({Ion2});
%
% atom=[atom, Ion1_atom, Ion2_atom];
%
% atom=atom_update_func(atom,MolID);
%
% All_resnames=sort(unique([atom.resname]));
% for i=1:length(All_resnames)
%     All_resnames(i)
% size(find(strcmp([atom.resname],All_resnames(i))),2)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script is used to translate all 'Resname' above/below certain ave_z
% ind_resname=find(strcmp([atom.resname],'TCH'))
% % ind_z=find([atom.ave_z]<5)
% % ind_resname=intersect(ind_resname,ind_z)
%
% x_shift=0.0;
% x_shift=num2cell([[atom(ind_resname).x]+x_shift]');
% [atom((ind_resname)).x]=deal(x_shift{:});
% z_shift=-.15;
% z_shift=num2cell([[atom(ind_resname).z]+z_shift]');
% [atom((ind_resname)).z]=deal(z_shift{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script is used to remove residues from a bilayer
%% in the simulation box
% keepnresidues=512;
% Box_dim=[150 150 Box_dim(3)]
% natomsinresidue=length(find([atom(:).molid]==1));
% ind_zlo=find([atom.med_z]<mean([atom.med_z]));
% ind_zhi=find([atom.med_z]>mean([atom.med_z]));
%
% ind_xhi=find([atom.med_x]<Box_dim(1));
% ind_yhi=find([atom.med_y]<Box_dim(2));
% ind_xy  = intersect(ind_xhi,ind_yhi);
% ind_zlo = intersect(ind_xy,ind_zlo);
% ind_zhi = intersect(ind_xy,ind_zhi);
%
% molid_zlo=unique([atom(ind_zlo).molid]);
% molid_zhi=unique([atom(ind_zhi).molid]);
%
% zlo_residues_left=length(ind_zlo)/natomsinresidue-keepnresidues/2;
% zhi_residues_left=length(ind_zhi)/natomsinresidue-keepnresidues/2;
%
% molid_zlo=molid_zlo(randperm(length(molid_zlo),zlo_residues_left));
% molid_zhi=molid_zhi(randperm(length(molid_zhi),zhi_residues_left));
%
% molid_rm=sort([molid_zlo molid_zhi]);
% rm_index=find(ismember([atom.molid],molid_rm));
% keep_ind=setdiff([ind_zlo ind_zhi],rm_index);
% atom=atom(keep_ind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% This script orders to atoms with respect to z
% sorted_coords=sortrows([[atom.index];[atom.molid];[atom.x];[atom.y];[atom.z];[atom.med_z]]',6);
% nAtoms=size(atom,2);
% Index=sorted_coords(:,1);
% MolID=sorted_coords(:,2);
% X_coord=sorted_coords(:,3);
% Y_coord=sorted_coords(:,4);
% Z_coord=sorted_coords(:,5);
%
% nmol=1;first_in=[1];last_in=[];
% for i=1:nAtoms;
%     if i > 1 && MolID(i) ~= MolID(i-1)
%         nmol=nmol+1;
%         atom(i).molid=nmol;
%         first_in(atom(i).molid,1)=i;last_in(atom(i).molid-1,1)=i-1;
%     elseif i > 1 && MolID(i) == MolID(i-1)
%         atom(i).molid=atom(i-1).molid;
%     elseif i == 1
%         atom(i).molid=1;
%     end
%
%     %atom(i).resname=Resname(i);
%     %atom(i).type=Atomtype(i);
%     %atom(i).fftype=Atomtype(i);
%     %atom(i).index=mod(i,100000);
%     %atom(i).neigh.type  = {};
%     atom(i).neigh.index  = zeros(6,1);
%     atom(i).neigh.dist  = zeros(6,1);
%     atom(i).bond.type  = zeros(6,1);
%     atom(i).bond.index  = zeros(6,1);
%     atom(i).x=X_coord(i);
%     atom(i).y=Y_coord(i);
%     atom(i).z=Z_coord(i);
% end
% last_in(atom(end).molid,1)=nAtoms;
%
% for i=1:max([atom(:).molid])
%     ind=find([atom.molid]==i);
%     [atom(ind).ave_x]=deal(mean([atom(ind).x]));
%     [atom(ind).ave_y]=deal(mean([atom(ind).y]));
%     [atom(ind).ave_z]=deal(mean([atom(ind).z]));
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script is used to remove residue
%% from the box, ie like and inverse 'replicate' function
%% Do not update atom...

% natomsinresidue=length(find([atom(:).molid]==1));
% nresidues=size(atom,2)/natomsinresidue;
% leftresidues=nresidues;
% newnresidues=64; % Edit manually
% i=1;nlo=0;nhi=0;rm=0;
% new_Box_dim=[1.1*Box_dim(1)/(nresidues/newnresidues)^.5 1.1*Box_dim(2)/(nresidues/newnresidues)^.5 Box_dim(3)];
% new_Box_dim=[75 75 95.75];
% while i<size(atom,2)
%             if mean([atom(i:(i-1)+natomsinresidue).z]) > Box_dim(3)/2 && nhi < newnresidues/2 && mean([atom(i:(i-1)+natomsinresidue).x]) < new_Box_dim(1) && mean([atom(i:(i-1)+natomsinresidue).y]) < new_Box_dim(2);
%                 nhi=nhi+1;
%                 i=i+natomsinresidue;
%             elseif mean([atom(i:(i-1)+natomsinresidue).z]) < Box_dim(3)/2 && nlo < newnresidues/2 && mean([atom(i:(i-1)+natomsinresidue).x]) < new_Box_dim(1) && mean([atom(i:(i-1)+natomsinresidue).y]) < new_Box_dim(2);
%                 nlo=nlo+1;
%                 i=i+natomsinresidue;
%             else
%                 ind=[atom.molid]==atom(i).molid;
%                 atom(ind)=[];
%                 rm=rm+1;
%             end
%             if size(atom,2)/natomsinresidue<leftresidues;
%             leftresidues=size(atom,2)/natomsinresidue
%             end
% end
%
% if rm == (nresidues-newnresidues); disp('Success!'); end
%
% for i=1:size(atom,2)
%     atom(i).molid=ceil(i/natomsinresidue);
%     atom(i).index=i;
% end
% nAtoms=size(atom,2);
%
% % x_shift=num2cell([[atom.x]-min([atom.x])]'); [atom(:).x]=deal(x_shift{:});
% % y_shift=num2cell([[atom.y]-min([atom.y])]'); [atom(:).y]=deal(y_shift{:});
% % z_shift=num2cell([[atom.z]-min([atom.z])]'); [atom(:).z]=deal(z_shift{:});
%
% min([atom.x])
% min([atom.y])
% min([atom.z])
% max([atom.x])
% max([atom.y])
% max([atom.z])
%
% Box_dim=[max([atom.x]) max([atom.y]) Box_dim(3)];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% And now update the atoms script...
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nAtoms=size(dataArray{:,5}(:),1);
% MolID = str2double((dataArray{:,1}));
% Resname = strtrim(dataArray{:,2});
% Atomtype = strtrim(dataArray{:,3});
% AtomID = dataArray{:,4};
% X_coord = dataArray{:,5};
% Y_coord = dataArray{:,6};
% Z_coord = dataArray{:,7};
% X_velo = dataArray{:,8};
% Y_velo = dataArray{:,9};
% Z_velo = dataArray{:,10};
% disp('.gro file imported')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nAtoms=size([atom.x],2);
atom=atom_update_func(atom);

% if length(MolID) < nAtoms;
%     MolID=[MolID;[MolID(end)+1:nAtoms]'];
% end
% % clear atom
% nmol=1;first_in=[1];last_in=[];
% for i=1:nAtoms;
%     if i > 1 && MolID(i) ~= MolID(i-1)
%         nmol=nmol+1;
%         atom(i).molid=nmol;
%         first_in(atom(i).molid,1)=i;last_in(atom(i).molid-1,1)=i-1;
%     elseif i > 1
%         atom(i).molid=atom(i-1).molid;
%     elseif i == 1
%         atom(i).molid=1;
%     end
%     %         atom(i).resname=Resname(i);
%     %         atom(i).type=Atomtype(i);
%     %         atom(i).fftype=Atomtype(i);
%     atom(i).index=mod(i,100000);
%     %         atom(i).neigh.type  = {};
%     %         atom(i).neigh.index  = zeros(6,1);
%     %         atom(i).neigh.dist  = zeros(6,1);
%     %         atom(i).bond.type  = zeros(6,1);
%     %         atom(i).bond.index  = zeros(6,1);
%     %         atom(i).x=X_coord(i);
%     %         atom(i).y=Y_coord(i);
%     %         atom(i).z=Z_coord(i);
%     %         atom(i).vx=X_velo(i);
%     %         atom(i).vy=Y_velo(i);
%     %         atom(i).vz=Z_velo(i);
% end
%
%%%%%% Write structure to file %%%%%%
if writepdb  > 0; write_atom_pdb(atom,Bond_index,Box_dim,filename_out); end % Not finished
if writemol2 > 0; write_atom_mol2(atom,Bond_index,Box_dim,filename_out); end % Not finished
if writegro  > 0; write_atom_gro(atom,Box_dim,filename_out); end
if writexyz  > 0; write_atom_xyz(atom,Box_dim,filename_out); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Plot the structure %%%%%%%%%%%
if plot == 1 && (writepdb+writemol2+writexyz+writegro)>0;
    if writemol2  > 0;  system(strcat(char({PATH2VMD}),char(strcat({' '},{strcat(filename_out,'.mol2')})))); % Use VMD to plot the simulation box with no cell size data
    elseif writepdb  > 0;  system(strcat(char({PATH2VMD}),char(strcat({' '},{strcat(filename_out,'.pdb')})))); % Use VMD to plot the simulation box
    elseif writegro  > 0;  system(strcat(char({PATH2VMD}),char(strcat({' '},{strcat(filename_out,'.gro')})))); % Use VMD to plot the simulation box
    else writexyz  > 0;  system(strcat(char({PATH2VMD}),char(strcat({' '},{strcat(filename_out,'.xyz')})))); end % Use VMD to plot the simulation box
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
