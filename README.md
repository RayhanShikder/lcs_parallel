# LCS Row Parallel (CPU)
This repository contains three parallel algorithms of the LCS approach in MPI, OpenMP, and hybrid MPI-OpenMP platforms.
For all the approaches, both the versions of the row-wise independent algorithm were developed.

The "data/" directory contains the required data sets of our study. Data sets of DNA sequences with different lenghts reside inside this directory.

There is a makefile inside each approaches. To run the codes, do the following steps:
1. Run the make (command: make) command from the corresponding directory
2. Submit the job written inside the "myjob" file (command: ./myjob).
3. An "output.txt" file will be generated with the calculated length of the LCS and the execution time of the code.

In the MPI approach, the number of processes can be set from the "myjob" file of the corresponding directory. The number of threads, for the OpenMP approach, can also be set from the "myjob" file of the corresponding directory. Likewise, both the number of threads, and the number of processes for the Hybrid approach can be set from the "myjob" file inside the corresponding directories.

To use this tool in your own machine, please follow the instructions from https://github.com/RayhanShikder/lcs_parallel/wiki


