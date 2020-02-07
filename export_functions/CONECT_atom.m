%% CONECT_atom.m
% * This function prints the CONECT records sometimes used in pdb files
%
%% Version
% 2.07
%
%% Contact
% Please report bugs to michael.holmboe@umu.se
%
%% Examples
% # CONECT_atom(atom,Box_dim,1.25,2.25) 
%
function CONECT_atom(atom,Box_dim,maxrshort,maxrlong) 

nAtoms=size(atom,2);
atom=bond_angle_atom(atom,Box_dim,maxrshort,maxrlong);
B=[Bond_index(:,1:2); Bond_index(:,2) Bond_index(:,1)];
b1=sortrows(B);
filename='conect';
fid = fopen(strcat(filename,'.pdb'), 'wt');
for i=1:max(b1(:,1))
    ind=find(b1(:,1)==i);
    b2=b1(ind,2);
    fprintf(fid,'CONECT%5i%5i%5i%5i%5i%5i%5i',[i;b2]);
    fprintf(fid,'\r\n');
end

fprintf(fid,'MASTER    %5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i%5i\r\n',[0    0    0    0    0    0    0    0 nAtoms    0 i    0]);
fprintf(fid,'END');

fclose(fid);

% % % COLUMNS         DATA TYPE        FIELD           DEFINITION
% % % ---------------------------------------------------------------------------------
% % %  1 -  6         Record name      "CONECT"
% % %  7 - 11         Integer          serial          Atom serial number
% % % 12 - 16         Integer          serial          Serial number of bonded atom
% % % 17 - 21         Integer          serial          Serial number of bonded atom
% % % 22 - 26         Integer          serial          Serial number of bonded atom
% % % 27 - 31         Integer          serial          Serial number of bonded atom
%%%%% We do not use these below...
% % % 32 - 36         Integer          serial          Serial number of hydrogen bonded atom
% % % 37 - 41         Integer          serial          Serial number of hydrogen bonded atom
% % % 42 - 46         Integer          serial          Serial number of salt bridged atom
% % % 47 - 51         Integer          serial          Serial number of hydrogen bonded atom
% % % 52 - 56         Integer          serial          Serial number of hydrogen bonded atom
% % % 57 - 61         Integer          serial          Serial number of salt bridged atom

% % % COLUMNS         DATA TYPE     FIELD          DEFINITION
% % % ----------------------------------------------------------------------------------
% % %  1 -  6         Record name   "MASTER"
% % % 11 - 15         Integer       numRemark      Number of REMARK records
% % % 16 - 20         Integer       "0"
% % % 21 - 25         Integer       numHet         Number of HET records
% % % 26 - 30         Integer       numHelix       Number of HELIX records
% % % 31 - 35         Integer       numSheet       Number of SHEET records
% % % 36 - 40         Integer       numTurn        deprecated
% % % 41 - 45         Integer       numSite        Number of SITE records
% % % 46 - 50         Integer       numXform       Number of coordinate transformation
% % %                                              records  (ORIGX+SCALE+MTRIX)
% % % 51 - 55         Integer       numCoord       Number of atomic coordinate records
% % %                                              records (ATOM+HETATM)
% % % 56 - 60         Integer       numTer         Number of TER records
% % % 61 - 65         Integer       numConect      Number of CONECT records
% % % 66 - 70         Integer       numSeq         Number of SEQRES records
% % % MASTER        0    0    0    0    0    0    0    0 1808    0 1808    0
