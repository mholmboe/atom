%% This function is used in the by the autofit_ff scheme
function ff = write_ff(ff,Atom_labels,parameters,varargin)

nOpt=size(Atom_labels,2);
nParam=size(parameters,1);

if nargin == 3
    n=1;
    for i=1:numel(Atom_labels)
        ind=find(strcmp([ff.type],Atom_labels(i)));
        [ff(ind).sigma_nm]=parameters(n,1);
        [ff(ind).e_kJmol]=parameters(n+1,1);
        n=n+2;
    end
    if nParam>2*nOpt
        for i=1:numel(Atom_labels)
            ind=find(strcmp([ff.type],Atom_labels(i)));
            try
                [ff(ind).charge]=parameters(n,1);
            catch
            end
            n=n+1;
        end
    end
else
    VAR=varargin{1};
    n=1;
    for i=1:numel(Atom_labels)
        ind=find(strcmp([ff.type],Atom_labels(i)));
        [ff(ind).sigma_nm]=strcat(VAR,num2str(n));
        [ff(ind).e_kJmol]=strcat(VAR,num2str(n+1));
        if (n+1)<nParam
            n=n+2;
        end
    end
    if nParam>2*nOpt
        for i=1:numel(Atom_labels)
            ind=find(strcmp([ff.type],Atom_labels(i)));
            try
                [ff(ind).charge]=strcat(VAR,num2str(n));
            catch
            end
            n=n+1;
        end
    end
end