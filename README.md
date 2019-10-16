# WaterSim
This program runs a molecular dynamics simulation of a cubic box of water molecules using an SPC model. 
It will account for LJ interactions and force-shifted coulombic interactions.
This program is equilibrated in a NVT simulation for 100ps.

gfortran equilWater.f90

The trajectory outputted by this equilibration is then continued by the analSim.f90. 
This will run an NVE simulation which will output information about velocity autocorrelation coefficients, 
radial distribution functions, and mean squared displacement.

gfortran analSim.f90

Both of this program will take constants from constants.f90 and subroutines from mdSubroutines.f90
