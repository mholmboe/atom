%% keep_resname.m
% * This function removes all but the resnames
% * Tested 15/04/2017
% * Please report bugs to michael.holmboe@umu.se


%% Examples
% * atom = keep_resname(atom,{'SOL'})

function atom = keep_resname(atom,resnames)

ind=[];
for i = 1:length(resnames)
    resnames(i)
    ind=[ind find(strcmp([atom.resname],resnames(i)))];
end
atom=atom(ind);
atom=update_atom(atom);
