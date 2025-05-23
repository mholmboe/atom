%% molecule_atom.m
% * This function will set the molecule ID (MolId), residue name (Resname)
% * of the atom struct, ie making the atom struct a single molecule, as
% * well as optionally setting the atomtype names to the resp. elements.
%
% * Molid is a integer number
% * Resname is a character string
% * Element is a boolean 1 (true) or 0 (false)
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = molecule_atom(atom) % Basic input arguments
% # atom = molecule_atom(atom,Molid) % Setting the Molid to an integer number
% # atom = molecule_atom(atom,Molid,Resname) % Also adding a residue name
% # atom = molecule_atom(atom,Molid,Resname,Element) % To set atomtype names to the resp. elements (1), or not (0)

function atom = molecule_atom(atom,varargin)

% Default values, may be overwritten below
Molid=1;
Resname='MOL';
Element=0;

if nargin==2
    Molid=varargin{1};
elseif nargin==3
    Molid=varargin{1};
    Resname=varargin{2};
elseif nargin==4
    Molid=varargin{1};
    Resname=varargin{2};
    Element=varargin{3};
end

if numel(Molid)==1
    [atom.molid]=deal(Molid);
elseif numel(Molid) == size(atom,2)
    for i=1:size(atom,2)
        atom(i).molid=Molid(i);
    end

end

if iscell(Resname)
    Resname=char(Resname);
end
if numel(Resname)>0
    [atom.resname]=deal({Resname});
end

if Element>0
    atom=element_atom(atom);
end
