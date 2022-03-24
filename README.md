# Diablo 3.0

This version is based on the original version of J. R. Taylor (cf. DIABLO2.0 ReadMe)
It includes updated FFTW and pHDF5 libraries, and uses the new MPI directives for shared-distributed allocated memory.


## Typical Run Workflow:
1. Create a new run directory
2. Copy the grid* and input* from an example in /setups and edit them as desired
3. Run /proc/init/create_grid.m MATLAB script to generate a grid.h5 file
4. Modify the subroutines "create_flow" and "create_th" in /for/flow.f90 (or use one of the pre-defined IC_Types) to define the desired ICs
5. Run >> /for/setup_run ${NprocVertical} ${NprocHorizontal}    --  This will create an executable called “diablo” in your run directory with the prescribed processor domain decomposition.
7. Run the code from inside the run directory.  For example, if using MPI, you could do “mpirun -np N ./diablo >output.dat &” which will run a job in the background and write the output to the file output.dat
8. Some post processing routines can be found in the /proc directory along with descriptions on their use.



## Environment Notes
- Requires FFTW v3.3+
- Requires parallel HDF5 v1.12+


