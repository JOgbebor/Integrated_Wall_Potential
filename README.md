# Integrated_Wall_Potential
Fit lennard-jones interaction parameters for LAMMPS wall potentials such that 
the walls mimic a graphitic slit pore of prescribed width (W).

The goal is to simulate a confined fluid in LAMMPS, but to avoid including the 
wall atoms explicitly, so that the overall properties of the system are those 
of the confined fluid, not the fluid and walls.
