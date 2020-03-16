# LAMMPS-utils
Collection of codes for modifying, improving, and analyzing LAMMPS input and output files. LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) is one of the popular molecular dynamics simulation codes.

* __addshell.c__: This code modifies a LAMMPS data file to prepare it for the adiabatic shell model by Mitchell and Fincham. Refer to the comments inside the code for more details. Please cite the following publication if you make a use of this code:

  *Khalkhali M, Ma X, Zhang H, Liu Q, “Bulk and surface properties of gypsum: A comparison between classical force fields and dispersion-corrected DFT calculations”, Computational Materials Science, 164, 8-16 (2019)* ([https://doi.org/10.1016/j.commatsci.2019.03.045](https://doi.org/10.1016/j.commatsci.2019.03.045))
  
* __MKConventor.c__: This code converts atomic configuration files written in different data formats to others. Data files supported in the current version are:
  * LAMMPS data files
  * DL_POLY CONFIG file
  * Materials Studio XSD files
  * XYZ files
