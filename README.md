# An OpenMP based tool for finding LCS of DNA sequence data
This repository contains three parallel algorithms of the LCS approach in MPI, OpenMP, and hybrid MPI-OpenMP platforms.
For all the approaches, both the versions of the row-wise independent algorithm were developed.
Finally, for general use of the tool, the OpenMP based implementation has been made user friendly and provided in the **Tool** directory.

The "data/" directory contains the required data sets of our study. Data sets of DNA sequences with different lenghts reside inside this directory.
### Experimental Use

#### Prerequisites
OpenMP and OpenMPI need to be installed in the machine before running the experiments. The versions of **gcc**, **OpenMPI**, and **OpenMP** should be: gcc 4.8.5 or later, OpenMPI version 1.10.7 or later, OpenMP version 3.1 or later.
### How To?
There is a makefile inside each approaches. To reproduce the results in an experimental setup in a cluster computer, do the following steps:
1. Run the `make` command from the corresponding directory.
2. Change the input file name in the `myjob` file (last line).
3. Submit the job written inside the "myjob" file (run `./myjob` from terminal).
4. An "output.txt" file will be generated with the calculated length of the LCS and the execution time of the code.

In the MPI approach, the number of processes can be set from the `myjob` file of the corresponding directory. The number of threads, for the OpenMP approach, can also be set from the `myjob` file of the corresponding directory. Likewise, both the number of threads, and the number of processes for the Hybrid approach can be set from the `myjob` file inside the corresponding directories.



### General Use
To use this tool in your own machine, please follow the instructions from https://github.com/RayhanShikder/lcs_parallel/wiki


