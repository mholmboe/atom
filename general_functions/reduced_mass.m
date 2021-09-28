%% reduced_mass.m
% * This function calculates the reduced mass.
%
%% Version
% 2.10
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% * reduced_mass = reduced_mass(Atom_label1)
% * reduced_mass = reduced_mass(Atom_label1,Atom_label2)

function reduced_mass = reduced_mass(Atom_label1,varargin)

if nargin == 1
    Atom_label2=Atom_label1;
else
    Atom_label2=varargin{1};
end

mass1=mass_atom(Atom_label1);
mass2=mass_atom(Atom_label2);

reduced_mass=(mass1*mass2)/(mass1+mass2);

