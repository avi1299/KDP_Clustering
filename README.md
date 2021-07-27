# KDP Clustering

To enable multithreading in certain functions, execute the below command on your terminal. Rerun this everytime you start a new terminal

`$ export OMP_NUM_THREADS = <NUMBER OF THREADS>`

For consumer CPUs, you can give upto half the threads in your system, after which performance may drop.

To compile the program , navigate to teh `bin` folder and then use the makefile as:

`$ make test`

This creates the binary test.

To see all the options, run,

`$ ./test -h`

To check if the binary works, run the binary on the ExampleData files or effectively check if the following commands work:

`$ make run`
`$ make xtc_1001_silent`

To clean the generated compilation files and binaries, run,

`$ make clean`
