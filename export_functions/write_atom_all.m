%% write_atom_all.m
% * This function tries to write various files for you. Works best for
% * systems designed for Clayff...
%
%% Version
% 2.11
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # write_atom_all(atom,Box_dim,filename) % Basic input arguments
% # write_atom_all(atom,Box_dim,filename,1.25,2.25) % rmaxshort and rmaxlong
% # write_atom_all(atom,Box_dim,filename,1.25,2.25,'clayff','spc/e') % rmaxshort and rmaxlong, ff and water model
%
function write_atom_all(atom,Box_dim,filename,varargin)


if nargin>3
    maxrshort=cell2mat(varargin(1))
    maxrlong=cell2mat(varargin(2))
else
    maxrshort=1.25
    maxrlong=1.25
end

if nargin>5
    ffname=char(varargin{3});
    if nargin>6
        watermodel=char(varargin{4});
    else
        disp('Unknown watermodel, will try SPC/E')
    end
    if strncmpi(ffname,'clayff',5)
        try
            atom=clayff_atom(atom,Box_dim,ffname,watermodel)
        catch
            warning('Problem with assigning the clayff atomtypes');
            pause(3)
        end
        try
            atom = charge_atom(atom,Box_dim,'clayff',watermodel)
        catch
            warning('Problem with setting the charges');
            pause(3)
        end
        
    elseif strcmpi(ffname,'interface')
        try
            atom=interface_atom(atom,Box_dim,ffname,watermodel);
        catch
            warning('Problem with assigning the interface atomtypes');
            pause(3)
        end
        try
            atom = charge_atom(atom,Box_dim,'interface',watermodel)
        catch
            warning('Problem with setting the charges');
            pause(3)
        end
    elseif strcmpi(ffname,'interface15')
        try
            atom=interface15_atom(atom,Box_dim,ffname,watermodel);
        catch
            warning('Problem with assigning the interface15 atomtypes');
            pause(3)
        end
        try
            atom = charge_atom(atom,Box_dim,'interface15',watermodel)
        catch
            warning('Problem with setting the charges');
            pause(3)
        end
        
    else
        disp('Unknown forcefield, will try clayff')
        try
            atom=clayff_atom(atom,Box_dim,ffname,watermodel)
        catch
            warning('Problem with assigning the clayff atomtypes');
            pause(3)
        end
        try
            atom = charge_atom(atom,Box_dim,'clayff',watermodel)
        catch
            warning('Problem with setting the charges');
            pause(3)
        end
        
    end
end

%write_atom_pdb(atom,Box_dim,filename) % Without CONECT records
try
    write_atom_pdb(atom,Box_dim,filename,maxrshort,maxrlong) % With CONECT records
catch
    warning('Problem with writing a .pdb file');
    pause(3)
end

try
    write_atom_gro(atom,Box_dim,filename)
catch
    warning('Problem with writing a .gro file');
    pause(3)
end

try
    write_atom_xyz(atom,Box_dim,filename)
catch
    warning('Problem with writing a .xyz file');
    pause(3)
end

try
    write_atom_cif(atom,Box_dim,filename)
catch
    warning('Problem with writing a .cif file');
    pause(3)
end

if nargin>5
    try
        write_atom_mol2(atom,Box_dim,filename,maxrshort,maxrlong,ffname,watermodel)
    catch
        warning('Problem with writing a .mol2 file');
        pause(3)
    end
    
    try
        write_atom_pqr(atom,Box_dim,filename,maxrshort,maxrlong,ffname,watermodel)
    catch
        warning('Problem with writing a .pqr file');
        pause(3)
    end
    
    try
        write_atom_psf(atom,Box_dim,filename,maxrshort,maxrlong,ffname,watermodel)
    catch
        warning('Problem with writing a .psf file');
        pause(3)
    end
    
    try
        write_atom_lmp(atom,Box_dim,filename,maxrshort,maxrlong,ffname,watermodel)
    catch
        warning('Problem with writing a lammps topology file');
        pause(3)
    end
    
    try
        write_atom_itp(atom,Box_dim,filename,maxrshort,maxrlong,ffname,watermodel)
    catch
        warning('Problem with writing a gromacs .itp file');
        pause(3)
    end
    
end

assignin('caller','out_atom',atom);
