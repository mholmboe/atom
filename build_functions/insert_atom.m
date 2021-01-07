%% insert_atom.m
% * This function inserts a whole molecule from a structure file or atom_in
% into a region defined by <limits> with a atom (molecule) structure
% * rotate can be a string like 'random', {'random'}, or be used to set
% some angles like [60 90 60]. varargin can be used to assure that one
% atom type is at least some distance above (in z) some other atom type.
%
%% Similar
% fuse_atom
% protonate_atom
% create_atom.m
%
%% Version
% 2.082
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * atom = insert_atom(atom,limits,'rotate',r,maxsol) % Basic input arguments
% * atom = insert_atom(atom,limits,[10 20 30],r,maxsol,solute_atom)
% * atom = insert_atom(atom,limits,'rotate',r,maxsol,solute_atom,{'C1' 'N1'},0.3)
% * atom = insert_atom(atom,limits,'rotate',r,maxsol,solute_atom,[1 4],0.3)

function atom = insert_atom(atom_in,limits,rotate,r,nmax,varargin)

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

if nargin>5
    solute_atom=varargin{1};
else
    solute_atom=[];
end

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

while (size(all_atom,2) < nmax*nAtoms_in) || n < 1000
    
    %     size(all_atom,2)
    %     nmax*nAtoms_in
    %     pause
    %     size(System)
    %     size(rotate_atom(atom_in,Box_dim,90*rand(1),90*rand(1),90*rand(1)))
    %     pause
    
    if ~mod(n,10) % Test of stride
        disp('% succesful attempts')
        (size(System,2)/nAtoms_in)/n*100
        if (size(System,2)/nAtoms_in)/n*100 < 1 %
            n=123456788;
        end
    end % end of test of stride
    
    temp_atom=rotate_atom(atom_in,limits(4:6),rotate);
    temp_atom=place_atom(temp_atom,[Lx*rand(1)+limits(1) Ly*rand(1)+limits(2) Lz*rand(1)+limits(3)]); % was translate_atom
    temp_atom=wrap_atom(temp_atom,limits(4:6));
    temp_atom=slice_atom(temp_atom,limits,0);
    if nargin>6
        difference=0;
        Atom_labels=varargin(2);
        type1=Atom_labels{1}(1)
        type2=Atom_labels{1}(2)
        if nargin>7
            difference=cell2mat(varargin(3));
        end
        if find(strcmp(type1,[temp_atom.type]))>0
            if median([temp_atom(strncmpi(type1,[temp_atom.type])).z],1) < (difference + median([temp_atom(strncmp1(type2,[temp_atom.type],1)).z]))
                temp_atom=[];
            else
                disp('Adding a molecule!!!')
            end
        else
            if [temp_atom(type1).z] < (difference + [temp_atom(type2).z])
                temp_atom=[];
            else
                disp('Adding a molecule!!!')
            end
        end
    end
    if size(temp_atom,2)>0
        try temp_atom = rmfield(temp_atom,'Mw');end
        try temp_atom = rmfield(temp_atom,'element');end
        try temp_atom = rmfield(temp_atom,'COM_x');end
        try temp_atom = rmfield(temp_atom,'COM_y');end
        try temp_atom = rmfield(temp_atom,'COM_z');end
        try assignin('caller','temp_atom',temp_atom);end
        try assignin('caller','all_atom',all_atom);end
        try assignin('caller','System_insert_atom',System);end
        try temp_atom = merge_atom(System,limits(4:6),temp_atom,'molid','H',[r-.5 r]);end
    end
    if size(temp_atom,2)>0
        System=wrap_atom(update_atom({System temp_atom}),limits(4:6));
        if size(all_atom,2)>0
            all_atom=update_atom({all_atom temp_atom});
        else
            all_atom=temp_atom;
        end
    end
    
    disp('Attempts Inserted nMax')
    [n nmax size(all_atom,2)/nAtoms_in]
    
    n=n+1;
    
    if size(all_atom,2) > nmax*nAtoms_in-nAtoms_in
        break
    end
end

disp('nMOL after merge')
size(all_atom,2)/nAtoms_in

% Do we need unwrapping here?
% all_atom = unwrap_atom(all_atom,limits(4:6),'xyz');

atom=all_atom;
if size(all_atom,2)>0
    atom=update_atom(all_atom);
else
    atom=atom_in;
    atom(1:end)=[];
end

if nmax==0
    try
        atom(1:end)=[];
    catch
    end
end
assignin('caller','limits',limits);



