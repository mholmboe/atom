%% keep_resname.m
% * This function removes all but the passed resnames
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # keep_resname(atom,{'SOL'})

function atom = keep_resname(atom,resnames)

ind=[];
for i = 1:length(resnames)
    resnames(i)
    ind=[ind find(strcmp([atom.resname],resnames(i)))];
end
atom=atom(ind);
atom=update_atom(atom);
