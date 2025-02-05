%% This special function is used to write the ffnonboned.itp file with VAR1..N string used in the ff optimization protocol
function ff = write_ffnonbonded(ffname,varargin)

if nargin>1
    Opt_labels=varargin{1};
    parameters=varargin{2};
else
    Opt_labels=[];
    parameters=[];
end

if nargin>3
    FEP_ion=varargin{3};
else
    FEP_ion=[];
end

if nargin>4
    FEP_charge=varargin{4};
else
    FEP_charge=[];
end

filename_out = 'ffnonbonded.itp';
file_title = 'Gromacs .itp topology file written in MATLAB from parameters in ...ff.mat'; % Header in output file

if ~isstruct(ffname)
    load(strcat(ffname,'.mat'));
    ff=ff.ff;
else
    ff=ffname;
    ffname='Forcefield parameters';
end


wat=load(strcat('water_models.mat'));
watff=wat.ff;

AllAtom_labels=[ff.type];
RestAtom_labels=AllAtom_labels;
if numel(Opt_labels)>0
    RestAtom_labels=AllAtom_labels(~ismember(AllAtom_labels,Opt_labels));
end
WatAtom_labels=[watff.type];

fid = fopen(filename_out, 'wt');
fprintf(fid, '%s\r\n','; Created in Matlab by reading the all atomypes (from a ff.mat file) and their associated parameters');
fprintf(fid, '%s % s\r\n',';',file_title);
fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','[ atomtypes ]');
fprintf(fid, '%s\r\n','; name  number  mass  charge  ptype  sigma   epsilon ;');

nOpt=size(Opt_labels,2);
nParam=size(parameters,1);

if numel(Opt_labels)>0
    ff = write_ff(ff,Opt_labels,parameters,'VAR');

    for i=1:numel(Opt_labels)
        ind=find(strcmp([ff.type],Opt_labels(i)));

        if  nParam>2*nOpt
            atomtypes(i,:) = {char([ff(ind).type]),[ff(ind).atnum],[ff(ind).mass],...
                char([ff(ind).charge]),'A',char([ff(ind).sigma_nm]),char([ff(ind).e_kJmol])};
            fprintf(fid, '%-12s % 3i\t% 9.5f\t% 9s   % -3s\t% -8s\t% -8s\r\n', atomtypes{i,:});
        else
            atomtypes(i,:) = {char([ff(ind).type]),[ff(ind).atnum],[ff(ind).mass],...
                [ff(ind).charge],'A',char([ff(ind).sigma_nm]),char([ff(ind).e_kJmol])};
            fprintf(fid, '%-12s % 3i\t% 9.5f\t% 9.6f   % -3s\t% -8s\t% -8s\r\n', atomtypes{i,:});
        end
    end

end

if numel(FEP_ion)>0
    ind=find(strncmpi([ff.type],FEP_ion,2));
    if numel(ind)==0
        ind=find(strncmpi([ff.type],FEP_ion,1));
    end

    if numel(ind)>1 && sum(ismember([ff(ind).charge],FEP_charge))>0
        ind=ind(abs([ff(ind).charge])==FEP_charge); % Ag+, Ag2+...
    else
        ind=ind(1);
    end
    ind(2)=ind;

    FEP_Atom_labels={'Io' FEP_ion  'Na' 'Cl'};FEP_Atom_labels=unique(FEP_Atom_labels,'stable')
    RestAtom_labels=RestAtom_labels(~ismember(RestAtom_labels,FEP_Atom_labels));
    for a=3:numel(FEP_Atom_labels)
        ind(a)=find(strncmpi([ff.type],FEP_Atom_labels{a},2));
    end
    fprintf(fid, '%s\r\n','; Current FEP ion (always named Io). Remember to also set its charge/s correctly in ions.itp');
    for i=1:numel(FEP_Atom_labels)
        atomtypes(i,:) = {char(FEP_Atom_labels(i)),[ff(ind(i)).atnum],[ff(ind(i)).mass],...
            [ff(ind(i)).charge],'A',[ff(ind(i)).sigma_nm],[ff(ind(i)).e_kJmol]};
        fprintf(fid, '%-12s % 3i\t% -9.5f\t% 9.6f   % -3s\t% -8.6e\t% -8.6e\r\n', atomtypes{i,:});
    end
end

fprintf(fid, '%s\r\n',strcat('; ',ffname));
for i=1:numel(RestAtom_labels)
    ind=find(strcmp([ff.type],RestAtom_labels(i)));
    atomtypes(i,:) = {char([ff(ind).type]),[ff(ind).atnum],[ff(ind).mass],...
        [ff(ind).charge],'A',round2dec([ff(ind).sigma_nm],6),round2dec([ff(ind).e_kJmol],6)};
    fprintf(fid, '%-12s % 3i\t% -9.5f\t% 9.6f   % -3s\t% -8.6e\t% -8.6e\r\n', atomtypes{i,:});
end

for i=1:numel(WatAtom_labels)
    ind=find(strcmp([watff.type],WatAtom_labels(i)));
    atomtypes(i,:) = {char([watff(ind).type]),[watff(ind).atnum],[watff(ind).mass],...
        [watff(ind).charge],'A',round2dec([watff(ind).sigma_nm],6),round2dec([watff(ind).e_kJmol],6)};
    fprintf(fid, '%-12s % 3i\t% -9.5f\t% 9.6f   % -3s\t% -8.6e\t% -8.6e\r\n', atomtypes{i,:});
end

fprintf(fid, '\r\n');

fclose(fid);

disp('ffnonbonded.itp structure file written')
