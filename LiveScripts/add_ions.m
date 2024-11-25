%% How-to add ions to a simulation box
% Add ions randomly wherever there is space (this function is similar to |solvate_atom()| 
% and |insert_atom()|
% 
% 

Box_dim=[30 30 30];
Cation='Na'
nCation=20;
System=[];
%% atom = create_atom(type,resname,limits,scale,maxion,in_atom)
Ion1 = create_atom(Cation,Cation,[0 0 0 Box_dim(1) Box_dim(2) Box_dim(3)],nCation/2,2,System);
% vmd(Ion,Box_dim)
show_atom(Ion1,Box_dim)

%%

%% Add ions using the copy_atom function
%% atom = copy_atom(atom,origtype,newtype,newresname,trans_vec,num);
Ion2 = copy_type(MMT,'Mgo',Cation,Cation,[0 0 8],sum(ismember([MMT.type],'Mgo'))/2);

%%

Ion3 = create_grid_atom('Na',nCatIonblk/planes,'Cl',0,[Box_dim(1) Box_dim(2) 95+i*8],'xy');