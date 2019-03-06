# Virtual FCS Measurements using MD Simulation Data

This goal of this notebook is to generate a simulated FCS measurement using data from a GROMACS simulation.

The setup of the virtual system is a membrane with a circularly symmetric incident beam.

Data is read in from an .xtc or .trr file, using a .gro file to define the system topology. Data frames from the input file are iterated through, and for each frame, a detected intensity from each lipid is calculated. The intensity trace and autocorrelation functions for each lipid, and for the total at each frame, are plotted.

## Dependencies

All dependencies are available via pip

- ```scipy```
- ```mdtraj```
- ```matplotlib```
- ```numpy```
- ```multiprocessing``` - Not technically required, but you'll be kicking yourself if you don't do the analysis in parallel
