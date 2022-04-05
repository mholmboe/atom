clear all;
format compact;
%% Import_the hexagonal Pyrophyllite structure, cutted and sliced in Avogadro for instance
%%
filename='hex_n19356_15.5nm_pyro_clayff.gro';
oldresname='pyro'; % New resname
newresname='MMT'; % New resname
import_atom(filename);

%% Find the edge Al's and rename them to a dummy name - DAl
%%
ind_Oalhhh=find(ismember([atom.type],{'Oalh' 'Oalhh'}));
fb_atom = find_bonded_atom(atom,Box_dim,ind_Oalhhh,'Al');
[atom(type2_ind).type]=deal({'DAl'});
%% Performing the iso subst, but not at the edge close to the DAl atoms
%%
NumOctSubst=2*ceil(floor(2/3*size(atom,2)/40)/2);
atom1 = substitute_atom(atom,Box_dim,NumOctSubst,'Al','Mgo',5);
[atom1(type2_ind).type]=deal({'Al'}); % Rename DAl atoms back to Al
atom1=clayff211_atom(atom1,Box_dim,'clayff','spc');
% vmd(atom1,Box_dim)
%% Write the gromacs itp file
%%
[atom1.resname]=deal({newresname});
filename=strrep(filename,oldresname,newresname);
write_atom_gro(atom1,Box_dim,filename)
write_atom_itp(atom1,Box_dim,erase(filename,'.gro'),1.25,2.25,'clayff','spc/e')
-NumOctSubst