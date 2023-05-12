# Integrated_Wall_Potential
This Python script writes a sub-script which calculates the integrated potential of atom walls with respect to some fluid confined between them.
It also fits Lennard Jones parameters (12-6, 10-4, or 9-3) to replicate the wall. Useful for simulating confined fluids without explicitly simulating the walls themselves.
Use the fitted parameters in the LAMMPS command: "fix wall/lj..." to mimic walls. 
The master file can execute the file it creates (last line of the script), which will produce relevant graphs. Ensure the exec(open().read()) command has the correct directory as its argument.
