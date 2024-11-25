clear all;clc;

AllAll_Neighbours=[];
atom=import_atom('1xMMT.gro');
atom=clayff_atom(atom,Box_dim);
AllAll_Neighbours=[AllAll_Neighbours All_Neighbours];

atom=import_atom('1xMMT_1octq.gro');
atom=clayff_atom(atom,Box_dim);
AllAll_Neighbours=[AllAll_Neighbours ; All_Neighbours];

atom=import_atom('1xMMT_1tetraq.gro');
atom=clayff_atom(atom,Box_dim);
AllAll_Neighbours=[AllAll_Neighbours ; All_Neighbours];

atom=import_atom('1xMMT_2octq.gro');
atom=clayff_atom(atom,Box_dim);
AllAll_Neighbours=[AllAll_Neighbours ; All_Neighbours];

atom=import_atom('1xMMT_3octtetraq.gro');
atom=clayff_atom(atom,Box_dim);
AllAll_Neighbours=[AllAll_Neighbours ; All_Neighbours];

atom=import_atom('Vermiculite_0000019.pdb');
atom=clayff_atom(atom,Box_dim);
AllAll_Neighbours=[AllAll_Neighbours ; All_Neighbours];

chargevec=[1.782 1.782 1.562 1.562 1.884 0.4];

AllNew_Neighbours=[];
atom=import_atom('1xMMT.gro');
atom=clayff_atom(atom,Box_dim);
atom=charge_clayff_atom(atom,Box_dim,{'Al' 'Alt' 'Mgh' 'Mgo' 'Si' 'H'},chargevec)
atom=ff_atom(atom,Box_dim);
AllNew_Neighbours=[AllNew_Neighbours All_Neighbours];

atom=import_atom('1xMMT_1octq.gro');
atom=clayff_atom(atom,Box_dim);
atom=charge_clayff_atom(atom,Box_dim,{'Al' 'Alt' 'Mgh' 'Mgo' 'Si' 'H'},chargevec)
atom=ff_atom(atom,Box_dim);
AllNew_Neighbours=[AllNew_Neighbours ; All_Neighbours];

atom=import_atom('1xMMT_1tetraq.gro');
atom=clayff_atom(atom,Box_dim);
atom=charge_clayff_atom(atom,Box_dim,{'Al' 'Alt' 'Mgh' 'Mgo' 'Si' 'H'},chargevec)
atom=ff_atom(atom,Box_dim);
AllNew_Neighbours=[AllNew_Neighbours ; All_Neighbours];

atom=import_atom('1xMMT_2octq.gro');
atom=clayff_atom(atom,Box_dim);
atom=charge_clayff_atom(atom,Box_dim,{'Al' 'Alt' 'Mgh' 'Mgo' 'Si' 'H'},chargevec)
atom=ff_atom(atom,Box_dim);
AllNew_Neighbours=[AllNew_Neighbours ; All_Neighbours];

atom=import_atom('1xMMT_3octtetraq.gro');
atom=clayff_atom(atom,Box_dim);
atom=charge_clayff_atom(atom,Box_dim,{'Al' 'Alt' 'Mgh' 'Mgo' 'Si' 'H'},chargevec)
atom=ff_atom(atom,Box_dim);
AllNew_Neighbours=[AllNew_Neighbours ; All_Neighbours];

atom=import_atom('Vermiculite_0000019.pdb');
atom=clayff_atom(atom,Box_dim);
atom=charge_clayff_atom(atom,Box_dim,{'Al' 'Alt' 'Mgh' 'Mgo' 'Si' 'H'},chargevec)
atom=ff_atom(atom,Box_dim);
AllNew_Neighbours=[AllNew_Neighbours ; All_Neighbours];


