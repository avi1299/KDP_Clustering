# KDP Clustering

This tool is designed to be a swiss army knife for analyzing the clustering of entities in Molecular Dynamic Simulations. Currently, it works on the Potassium Dihydrogen Phosphate(KDP) system, but can be adapted to other systems by replacing a few files.

***
## Features
The tool has the following features:
1. Ability to read either **.PDB** files or **.XTC** and **.TOP** files
2. A .TOP file parser
3. Clustering of HPO molecules using the Depth-First-Search algorithm [Time complexity ~ O(n)]
4. Finding associated counterion K molecules associated with each cluster
5. Drawing out statistics on charge and size of each cluster
6. Calculating the coordination number of the cluster
7. Ability to switch on and switch off the periodic boundary condition
8. Ring analysis on the entire configuration to identify rings
[Time complexity ~ O(n<sup>3</sup>)]
9. Multithreaded operation in some functions
10. Multilevel verbose output
11. Output of HPO clusters and K ions according to the cluster sizes
12. Ability to differentiate between strong and weak bonds between molecules
13. Captures the directionality of bonds between molecules
13. Ability to specify stricter conditions for connection between two molecules and checking the difference in the number of strict and normal connections
14. Ability to output statistics in terms of percentage and numbers
15. Ability to specify whether the input is a solution or a supercell

These are the planned features:
1. Calculate Center of Mass, Radius of gyration of the clusters
2. Track cluster through the simulation and record when molecules join or leave the cluster
3. Find the distribution of charge and mass in the cluster via quadrupole and higher moments
4. Graph the quantities as they evolve with time
5. Report Potassium association to cluster wrt time
6. Define and implement a metric for layer formation in the configuration
7. Multithread across configurations

***
## Requirements
We require the GNU C++ compiler `g++` to compile our code into a binary. The included GROMACS library requires the BLAS and LAPACK libraries to run.
It is suggested you install the below mentioned libraries.

1. `g++ >= 9.3.0`
2. `libblas-dev >= 3.9.0`
3. `liblapack-dev >= 3.9.0`
***
## Setup

Clone the repo to your system and change the directory to the folder:

`$ git clone https://github.com/avi1299/KDP_Clustering.git`

`$ cd KDP_Clustering`

To enable multithreading in certain functions, execute the below command on your terminal. Rerun this every time you start a new terminal

`$ export OMP_NUM_THREADS = <NUMBER OF THREADS>`

For consumer CPUs, you can give up to half the threads in your system, after which performance may drop.

To compile the program, navigate to the `bin` folder 

`$ cd bin`

and then use the makefile as:

`$ make test`

This creates the binary test. To see all the options, run,

`$ ./test -h`

```
Usage: ./test [OPTIONS].. [ARGUMENTS]..
  -h            : Prints the help information and exits
  -f <file>     : Specify the input file
  -t <file>     : Specify the .top file in case of XTC input
  -o <file>     : Specify the output file
  -v <int>      : Specify the verbose level among {1,2,3,4}.
  -c            : Enables checking if stricter connectedness conditions affects number of connections in each configuration
  -s <int>      : Specify size of clusters to be outputted in PDB format
  -g            : Used with -s to include clusters having size greater than or equal to the argument for -s
  -p            : Prints the percentage of molecule belonging to a cluster rather than the number of molecules
  -m <int>      : Prints statisitics every argument number of configurations(default=1). Needs verbose flag to be disabled
  -n            : Uses only strong bonds to perform clustering
  -r            : Performs ring analysis and outputs rings instead
  -a            : Specifies that the input file is that of a supercell and therefore Periodic Boundary Condition will not applied
```

To check if the binary works, run the binary on the ExampleData files.
The makefile has commands for the same:

`$ make run`

`$ make xtc_1001_silent`

These commands will run the binary and will produce the output files
`test_out.pdb`, `cluster_max_size.dat`, `cluster_statistics.dat` and `ring_statistics.dat` if ring analysis is enabled. 

To clean the generated compilation files and binaries, run,

`$ make clean`

***




