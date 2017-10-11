%% insert_atom.m
% * This function inserts a molecule from a structure file or atom_in into 
% * a region defined by <limits> with a atom (molecule) structure
% * rotate can be a string like 'random', {'random'}, or be used to set
% * some angles like [60 90 60]. varargin can be used to assure that one
% * atom type is at least some distance above (in z) some other atom type.
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se

%% Examples
% * atom = insert_atom(atom,limits,rotate,r,maxsol,solute_atom)
% * atom = insert_atom(atom,limits,rotate,r,maxsol,solute_atom,{'C1'},{'C2'},0.3)

function atom = insert_atom(atom_in,limits,rotate,r,nmax,solute_atom,varargin)

if numel(limits)==1
    Lx=limits(1);
    Ly=limits(1);
    Lz=limits(1);
    limits(4)=limits(1);
    limits(5)=limits(1);
    limits(6)=limits(1);
    limits(1:3)=0;
elseif numel(limits)==3
    Lx=limits(1);
    Ly=limits(2);
    Lz=limits(3);
    limits(4)=limits(1);
    limits(5)=limits(2);
    limits(6)=limits(3);
    limits(1:3)=0;
elseif numel(limits)==6
    Lx=limits(4)-limits(1);
    Ly=limits(5)-limits(2);
    Lz=limits(6)-limits(3);
end

nAtoms_in=size(atom_in,2);

n=0;
if size(solute_atom,2) == 0
    all_atom=center_atom(atom_in,[Lx Ly Lz],'all','xyz');
    all_atom=translate_atom(all_atom,[limits(1) limits(2) limits(3)],'all');
    System=wrap_atom(all_atom,limits(4:6));
    n=1;
else
    System=wrap_atom(solute_atom,limits(4:6));
    all_atom=[];
end

while size(all_atom,2) < nmax*nAtoms_in-nAtoms_in || n < 1000
    %     size(System)
    %     size(rotate_atom(atom_in,Box_dim,90*rand(1),90*rand(1),90*rand(1)))
    %     pause
    n
    temp_atom=rotate_atom(atom_in,limits(4:6),rotate);
    temp_atom=translate_atom(temp_atom,[Lx*rand(1) Ly*rand(1) Lz*rand(1)],'all');
    temp_atom = wrap_atom(temp_atom,limits(4:6));
    temp_atom = slice_atom(temp_atom,limits,0);
    if nargin>6
        difference=0;
        Atom_labels=varargin(1);
        type1=Atom_labels{1}(1);
        type2=Atom_labels{1}(2);
        if nargin>7
            difference=cell2mat(varargin(2));
        end
        if median([temp_atom(strncmpi(type1,[temp_atom.type],1)).z]) < (difference + median([temp_atom(strncmpi(type2,[temp_atom.type],1)).z]));
            temp_atom=[];
        end
    end
    if size(temp_atom,2)>0
        temp_atom = rmfield(temp_atom,'Mw');
        temp_atom = rmfield(temp_atom,'element');
        temp_atom = rmfield(temp_atom,'COM_x');
        temp_atom = rmfield(temp_atom,'COM_y');
        temp_atom = rmfield(temp_atom,'COM_z');
        assignin('caller','temp_atom',temp_atom);
        assignin('caller','all_atom',all_atom);
        assignin('caller','Systemf',System);
        temp_atom = merge_atom(System,limits(4:6),temp_atom,'molid','H',[r-.5 r]);
    end
    if size(temp_atom,2)>0
        all_atom=update_atom({all_atom temp_atom});
        System=wrap_atom(update_atom({System temp_atom}),limits(4:6));
    end
    n=n+1;
    
    if size(all_atom,2) > nmax*nAtoms_in-nAtoms_in
        break
    end
end

disp('nMOL after merge')
size(all_atom,2)/nAtoms_in

all_atom = unwrap_atom(all_atom,limits(4:6),'xyz');

atom=all_atom;
if size(all_atom,2)>0
    atom=update_atom(all_atom);
end
assignin('caller','limits',limits);



