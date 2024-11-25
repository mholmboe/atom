%% import_mclf_C6dispersion.m
% * This function imports the C6 dispersion terms from the
% mclf code. Not tested alot..
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # atom = import_atom_ddec('XYZ_even_tempered_net_atomic_charges.xyz')
%
function All_CX = import_mclf(varargin)

if nargin>0
    filename=varargin{1};
else
    filename='MCLF_screened_C6_dispersion_coefficients.xyz';
end

All_CX=[];
try
    [atomC6,Box_dim] = import_mclf_C6(filename);
    assignin('caller','All_C6',All_C6);
    All_CX=[All_CX;cell2mat(All_C6(:,3)')];
catch
    disp('Could not import the C6 dispersion results!')
end

try
    atomC8 = import_mclf_C8;
    assignin('caller','All_C8',All_C8);
    All_CX=[All_CX;cell2mat(All_C8(:,3)')];
catch
    disp('Could not import the C8 dispersion results!')
end

try
    atomC10 = import_mclf_C10;
    assignin('caller','All_C10',All_C10);
    All_CX=[All_CX;cell2mat(All_C10(:,3)')];
catch
    disp('Could not import the C10 dispersion results!')
end

try
    [atomq,Box_dim] = import_ddec_charges;
    assignin('caller','All_Q',All_Q);
    All_CX=[All_CX;cell2mat(All_Q(:,2)')];
catch
    disp('Could not import the DDEC charges!')
end

if exist("atomq","var")
    atom_mclf=atomq;
    [atom_mclf.charge]=atomq.charge;
end

if exist("atomC6","var")
    [atom_mclf.C6_aubhr6]=atomC6.C6_aubhr6;
    [atom_mclf.C6_kJmolnm6]=atomC6.C6_kJmolnm6;
end
if exist("atomC8","var")
    [atom_mclf.C8_aubhr8]=atomC8.C8_aubhr8;
    [atom_mclf.C8_kJmolnm8]=atomC8.C8_kJmolnm8;
end
if exist("atomC10","var")
    [atom_mclf.C10_aubhr10]=atomC10.C10_aubhr10;
    [atom_mclf.C10_kJmolnm10]=atomC10.C10_kJmolnm10;
end

Box_dim_mclf=Box_dim;
Atom_labels=unique([atom_mclf.type]);
assignin('caller','Atom_labels',Atom_labels);
assignin('caller','atom_mclf',atom_mclf);
assignin('caller','Box_dim_mclf',Box_dim_mclf);