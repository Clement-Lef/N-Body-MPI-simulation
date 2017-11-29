Parallel implementation of the N-Body problem applied to the gravity force.
The algorithm used is Barnes-Hut and the forces computations are parallelized using the MPI library.

To compile the program :

- Run the CMakeList with cmake .
- Run the makefile with make


To run the program the following arguments are necessary :

- Number of bodies
- Number of cluster of bodies
- Time step
- Final time
- Threshold for force computation
- Number of time step between each output

These arguments must be given in the order given above, example :

mpirun -np 2 ./NBodyMPI 40 1 0.01 50 3 1


The output is found in the file Output.dat in the following format for 3 time steps :

Body1_position_x Body1_position_y Body2_position_x Body2_position_y
Body1_position_x Body1_position_y Body2_position_x Body2_position_y
Body1_position_x Body1_position_y Body2_position_x Body2_position_y


The MATLAB file DrawBodies.m can be run to plot the output.

