%% gmx_make_ndx.m
% * This function helps you print custom gromacs .ndx files
%
%% Version
% 2.0
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% * gmx_make_ndx(groupname,index_vector)

function gmx_make_ndx(groupname,ind)

%clear all; format compact;
%You could set the groupname here manually if you want
%groupname='FIX_MMT_2';
%ind=[1+2880:40:2880+1992];
%This is in case you want to create a group with multiple selections
%ind1=1:10:960;
%ind2=6:10:960;
%ind3=10:10:960;
%ind=sort([ind1 ind2 ind3]);

% Format the ind vector to have 15 entries per row
ext_ind=zeros(1,15*ceil(length(ind)/15));
ext_ind(1:length(ind))=ind;
ext_ind=reshape(ext_ind,15,[])';

% Print the index file
fid = fopen(strcat(groupname,'.ndx'), 'wt');
fprintf(fid, '%s %s %s\r\n','[',groupname,']');
for i = 1:size(ext_ind,1);
    row=ext_ind(i,:);
    row(row==0)=[];
    if max(row)<1000;
    fprintf(fid, '%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i%4i\r\n', row);
    elseif max(row)<10000;
    fprintf(fid, '%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i\r\n', row);
    elseif max(row)<100000;
    fprintf(fid, '%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i\r\n', row);
    end
end
fprintf(fid, '\r\n');
fclose(fid);