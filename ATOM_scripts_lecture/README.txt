This folder contains some basic presentations about the ATOM library and an introduction to Gromacs, as well as simple example systems of clay systems using a the general and tailored MINFF forcefields, a modifed CLAYFF (containing new atomtypes/names), the original CLAYFF (Cygan, 2004) as well as a custom implementation of the Heinz, 2005 version of the Interface FF as well as the later Interface FF package from Heinz, 2013, using structures from the model database in the Interface FF 1.5 package.

Regarding the simulation files, if you are new to Gromacs and do not know where to start inspect the files called job.sh as well as topol.top and go from there! The systems have been tested with Gromacs 2018. More recent versions may demand using a -DFLEXIBLE flag in the .mdp files

NOTE: I cannot guarantee that the forcefield parameters and settings are correct, you should check this yourself by reading the original papers! 

Good papers to start with...
Cygan, R. T., Liang, J. J., & Kalinichev, A. G. (2004). Molecular models of hydroxide, oxyhydroxide, and clay phases and the development of a general force field. Journal of Physical Chemistry B, 108(4), 1255–1266.

Heinz, H., Koerner, H., Anderson, K. L., Vaia, R. A., & Farmer, B. L. (2005). Force Field for Mica-Type Silicates and Dynamics of Octadecylammonium Chains Grafted to Montmorillonite. Chemistry of Materials, (5), 5658–5669.
Heinz, H., Lin, T. J., Kishore Mishra, R., & Emami, F. S. (2013). Thermodynamically consistent force fields for the assembly of inorganic, organic, and biological nanostructures: The INTERFACE force field. Langmuir, 29(6), 1754–1765. https://doi.org/10.1021/la3038846

