%% write_atom_sdf.m
% * This function writes an sdf file from the atom struct
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_sdf(atom,Box_dim,filename_out) % Basic input arguments
%
function write_atom_sdf(atom,Box_dim,filename_out,varargin)

if regexp(filename_out,'.sdf') ~= false
    filename_out = filename_out;
else
    filename_out = strcat(filename_out,'.sdf');
end

if numel(Box_dim)==1
    Box_dim(1)=Box_dim(1);
    Box_dim(2)=Box_dim(1);
    Box_dim(3)=Box_dim(1);
end

% for i=1:size(atom,2)
%     if strncmpi(atom(i).type,{'Si'},2);atom(i).element={'Si'};atom(i).formalcharge=4;
%     elseif strncmpi(atom(i).type,{'SY'},2);atom(i).element={'Si'};atom(i).formalcharge=4;
%     elseif strncmpi(atom(i).type,{'SC'},2);atom(i).element={'Si'};atom(i).formalcharge=4;
%     elseif strncmpi(atom(i).type,{'S'},1);atom(i).element={'S'};atom(i).formalcharge=2;
%     elseif strncmpi(atom(i).type,{'AC'},2);atom(i).element={'Al'};atom(i).formalcharge=3;
%     elseif strncmpi(atom(i).type,{'AY'},2);atom(i).element={'Al'};atom(i).formalcharge=3;
%     elseif strncmpi(atom(i).type,{'Al'},2);atom(i).element={'Al'};atom(i).formalcharge=3;
%     elseif strncmpi(atom(i).type,{'Mg'},2);atom(i).element={'Mg'};atom(i).formalcharge=2;
%     elseif strncmpi(atom(i).type,{'Fe'},2);atom(i).element={'Fe'};atom(i).formalcharge=3;
%     elseif strncmpi(atom(i).type,{'O'},1);atom(i).element={'O'};atom(i).formalcharge=-2;
%     elseif strncmpi(atom(i).type,{'H'},1);atom(i).element={'H'};atom(i).formalcharge=1;
%     elseif strcmpi(atom(i).type,{'N'});atom(i).element={'N'};atom(i).formalcharge=0;
%     elseif strncmpi(atom(i).type,{'Ni'},2);atom(i).element={'Ni'};atom(i).formalcharge=2;
%     elseif strncmpi(atom(i).type,{'Li'},2);atom(i).element={'Li'};atom(i).formalcharge=1;
%     elseif strncmpi(atom(i).type,{'Na'},2);atom(i).element={'Na'};atom(i).formalcharge=1;
%     elseif strncmpi(atom(i).type,{'K'},1);atom(i).element={'K'};atom(i).formalcharge=1;
%     elseif strncmpi(atom(i).type,{'Cs'},2);atom(i).element={'Cs'};atom(i).formalcharge=1;
%     elseif strncmpi(atom(i).type,{'Co'},2);atom(i).element={'Co'};atom(i).formalcharge=2;
%     elseif strncmpi(atom(i).type,{'Cr'},2);atom(i).element={'Cr'};atom(i).formalcharge=3;
%     elseif strncmpi(atom(i).type,{'Cu'},2);atom(i).element={'Cu'};atom(i).formalcharge=2;
%     elseif strncmpi(atom(i).type,{'F'},1);atom(i).element={'F'};atom(i).formalcharge=-1;
%     elseif strncmpi(atom(i).type,{'Cl'},2);atom(i).element={'Cl'};atom(i).formalcharge=-1;
%     elseif strncmpi(atom(i).type,{'Br'},2);atom(i).element={'Br'};atom(i).formalcharge=-1;
%     elseif strncmpi(atom(i).type,{'I'},1);atom(i).element={'I'};atom(i).formalcharge=-1;
%     elseif strncmpi(atom(i).type,{'Ca'},2);atom(i).element={'Ca'};atom(i).formalcharge=2;
%     elseif strncmpi(atom(i).type,{'Sr'},2);atom(i).element={'Sr'};atom(i).formalcharge=2;
%     elseif strncmpi(atom(i).type,{'Ba'},2);atom(i).element={'Ba'};atom(i).formalcharge=2;
%     elseif strncmpi(atom(i).type,{'C'},1);atom(i).element={'C'};atom(i).formalcharge=0;
%     elseif strncmpi(atom(i).type,{'P'},1);atom(i).element={'P'};atom(i).formalcharge=0;
%     else
%         [atom(i).element]=atom(i).type;atom(i).formalcharge=0;
%     end
% end


if nargin>3
    long_r=cell2mat(varargin(1));
else
    long_r=2.25;
end
atom=bond_atom(atom,Box_dim,long_r,0.6);
nAtoms=size(atom,2);
nBonds=size(Bond_index,1);
assignin('caller','Bond_index',Bond_index);
% assignin('caller','atom_Bond_index',atom);

Atom_section=cell(nAtoms,10);
fid = fopen(filename_out, 'wt');

fprintf(fid, '%s\r\n',filename_out);
fprintf(fid, '%s\r\n','The atom MATLAB library, version 2.07');
fprintf(fid, '%s\r\n','Simple .sdf coordinate file with bonds');

%%%   9  8  0     0  0  0  0  0  0999 V2000
%%% 123456789012345678901234567890123456789
%%% The Counts line
Counts={nAtoms, nBonds, 0, 0, 0, 0, 0, 0, 0, 999, 'V2000'};
fprintf(fid, '%3i%3i%3i   %3i%3i%3i%3i%3i%3i%3i%6s\r\n',Counts{1:11});


%%% The Atoms block
%%% After the counts line comes the atoms block. For each atom mentioned in
%%% the first field of counts, include a line like so:
%%%    0.5369    0.9749    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
%%% The first three fields, 10 characters long each, describe the atom's
%%% position in the X, Y, and Z dimensions. After that there is a space, and
%%% three characters for an atomic symbol (O for oxygen, in this instance).
for i = 1:nAtoms
    if size(atom(i).type{1},2) > 3
        disp('Hey, this atom type name is actually too long for sdf')
        disp('chopping it down to 3 characters')
        [atom(i).index atom(i).type]
        atom(i).type=atom(i).type{1}(1:3);
    end
    Atomsblock(1:16)={atom(i).x, atom(i).y, atom(i).z, char(atom(i).type), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    fprintf(fid, '%10.4f%10.4f%10.4f %-3s%2i%3i%3i%3i%3i%3i%3i%3i%3i%3i%3i%3i\r\n',Atomsblock{1:16});
end

%%% The Bonds block
for i = 1:nBonds
    Bondsblock(1:7)={Bond_index(i,1), Bond_index(i,2), 1, 0, 0, 0, 0};
    fprintf(fid, '%i  %i  %i  %i  %3i%3i%3i\r\n',Bondsblock{1:7});
end

fprintf(fid,'M  END');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

fprintf(fid,'$$$$');

fprintf(fid, '\r\n');
fprintf(fid, '\r\n');

fclose(fid);

disp('.sdf structure file written')

