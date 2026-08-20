This folder contains forcefield parameters for general (G)MINFF, as well as tailored (T)MINFF parameters for specific mineral (see list below).

**Lennard-Jones parameters for the general MINFF**

Lennard-Jones parameters for all metal sites optimized for over 30 mineral synchronously, 
using the Lennard-Jones parameters for the oxygen atomtypes taken from OPC3 water model fixed.  

* Filename ffnonbonded_gminff.itp:
	- General MINFF for k=0, 250, 500, 1500 kJ/mol/rad2
	- CLAYFF
	- Several ion-pair potentials and water models


**Lennard-Jones parameters for the tailored MINFF**

Note that the Lennard-Jones parameters inn step 1 for each mineral were first optimized for only the metal sites only.
In a second step the oxygen Lennard-Jones parameters were optimized, keeping the metal parameters from step 1 fixed. 
These new oxygen parameters optimized in step 2 are commented by a ; and not necessarily better than the parameters optmized in step 1.

* Filename ffnonbonded_tminff_k0|k250|k500|1500.itp
	- tailored MINFF for k=0, 250, 500, 1500 kJ/mol/rad2, sorted by mineral

* Filename ffnonbonded_tminff_all_k_sorted_by_mineral.itp
	- tailored MINFF for k=0, 250, 500, 1500 kJ/mol/rad2, sorted by mineral
	
Mineral list (# 13 and 22 missing)
1	Kaolinite
2	Pyrophyllite
3	Talc
4	Forsterite
5	Brucite
6	Corundum
7	Quartz
8	Gibbsite
9	Li2O
10	Coesite
11	Cristobalite
12	Maghemite
14	Akdalaite
15	Boehmite
16	Diaspore
17	Periclase
18	Goethite
19	Hematite
20	Lepidocrocite
21	Wustite
23	CaF2
24	CaO
25	Portlandite
26	Nontronite
27	Montmorillonite
28	Dickite
29	Hectorite-F
30	Hectorite-H
31	Nacrite
32	Imogolite
33	Anatase
34	Rutile
35	cis_Oct_Fe2_cis
36	cis_Oct_Fe2_trans
37	cis_Oct_Mg2cis_Fe3cis
38	cis_Oct_Mg2cis_Fe3trans
39	cis_Oct_Mg2trans_Fe3cis
40	cis_Oct_Mg2trans_Fe3trans
41	cis_Tet_Fe3
42	trans_Oct_Fe2_cis
43	trans_Oct_Mg2cis_Fe3cis
44	trans_Tet_Fe3
45	Muscovite
