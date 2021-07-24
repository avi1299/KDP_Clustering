# KDP Clustering

To enable multithreading in certain functions, execute the below command on your terminal. Rerun this everytime you start a new terminal

`$ export OMP_NUM_THREADS = <NUMBER OF THREADS>`

For consumer CPUs, you can give upto half the threads in your system, after which performance may drop.

To compile the program use the makefile as:

`$ make test`

This creates the binary test.

To see all the options, run,

`$ ./test -h`

To clean the generated compilation files and binaries, run,

`$ make clean`
