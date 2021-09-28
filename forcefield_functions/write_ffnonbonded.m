%% This special function is used to write the ffnonboned.itp file with VAR1..N string used in the ff optimization protocoll 
function ff = write_ffnonbonded(Opt_labels,parameters,varargin)

if nargin>2
   ffname=varargin{1}; 
else
    ffname='ff';
end

filename_out = 'ffnonbonded.itp';
file_title = 'Gromacs .itp topology file written in MATLAB from parameters in ...ff.mat'; % Header in output file

ff=load(strcat(ffname,'.mat'));
ff=ff.ff;

wat=load(strcat('water_models.mat'));
watff=wat.ff;

AllAtom_labels=[ff.type];
RestAtom_labels=AllAtom_labels(~ismember(AllAtom_labels,Opt_labels));

WatAtom_labels=[watff.type];

fid = fopen(filename_out, 'wt');
fprintf(fid, '%s\n','; Created in Matlab by reading the all atomypes (from a ff.mat file) and their associated parameters');
fprintf(fid, '%s % s\r\n',';',file_title);
fprintf(fid, '\r\n');
fprintf(fid, '%s\r\n','[ atomtypes ]');
fprintf(fid, '%s\r\n','; name  number  mass  charge  ptype  sigma   epsilon ;');

nOpt=size(Opt_labels,2);
nParam=size(parameters,1);
ff = write_ff(ff,Opt_labels,parameters,'VAR');

for i=1:numel(Opt_labels)
    ind=find(strcmp([ff.type],Opt_labels(i)));
    atomtypes(i,:) = {char([ff(ind).type]),[ff(ind).atnum],[ff(ind).mass],...
        char([ff(ind).charge]),'A',char([ff(ind).sigma_nm]),char([ff(ind).e_kJmol])};
    if  nParam>2*nOpt
        fprintf(fid, '%-12s % 3i\t% 9.6f\t% 9s % 9s\t% 9s\t% 9s\r\n', atomtypes{i,:});
    else
        fprintf(fid, '%-12s % 3i\t% 9.6f\t% 9.6f % 9s\t% 9s\t% 9s\r\n', atomtypes{i,:});
    end
end

for i=1:numel(RestAtom_labels)
    ind=find(strcmp([ff.type],RestAtom_labels(i)));
    atomtypes(i,:) = {char([ff(ind).type]),[ff(ind).atnum],[ff(ind).mass],...
        [ff(ind).charge],'A',[ff(ind).sigma_nm],[ff(ind).e_kJmol]};
    fprintf(fid, '%-12s % 3i\t% 9.6f\t% 9.6f % 9s\t% 9.6g\t% 9.6g\r\n', atomtypes{i,:});
end

for i=1:numel(WatAtom_labels)
    ind=find(strcmp([watff.type],WatAtom_labels(i)));
    atomtypes(i,:) = {char([watff(ind).type]),[watff(ind).atnum],[watff(ind).mass],...
        [watff(ind).charge],'A',[watff(ind).sigma_nm],[watff(ind).e_kJmol]};
    fprintf(fid, '%-12s % 3i\t% 9.6f\t% 9.6f % 9s\t% 9.6g\t% 9.6g\r\n', atomtypes{i,:});
end

fprintf(fid, '\r\n');

fclose(fid);

disp('ffnonbonded.itp structure file written')
