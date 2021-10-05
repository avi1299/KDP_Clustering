IDIR =./include
SDIR =./src
ODIR=./obj
LDIR =./lib
CC=g++
CFLAGS=-I$(IDIR) -DGMX_DOUBLE=0 -fopenmp -mcmodel=medium #-fno-stack-protector#-pg -g
MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
MKFILE_DIR := $(dir $(MKFILE_PATH))
LD_LIBRARY_PATH=$(MKFILE_DIR)lib

LIBS=-lgromacs  -lm

_DEPS = DFS_clustering.h HPO.h input.h stack.h data_manip.h output.h graph.h constants.h charge.h bounding_box.h ring.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = HPO.o stack.o DFS_clustering.o input.o data_manip.o output.o graph.o charge.o bounding_box.o ring.o test.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))



$(ODIR)/%.o:$(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -L$(LDIR) $(LIBS) -Wl,-rpath=$(LD_LIBRARY_PATH)

test:$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -L$(LDIR) $(LIBS) -Wl,-rpath=$(LD_LIBRARY_PATH)

run: test
	./test -f ../ExampleData/M1A_S2_16.pdb -o ../ExampleData/test_out.pdb -v1

small: test
	./test -f ../ExampleData/M1A_S1_Avi.xtc -t ../ExampleData/M1A_S1.top -o ../ExampleData/test_out.pdb -m 100

#xtc_1: test
#	./test -f ../ExampleData/M2A_S3_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -v2

#xtc_1_silent: test
#	./test -f ../ExampleData/M2A_S3_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -m 10

xtc_1001: test
	./test -f ../ExampleData/M2A_S3_1001_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -v2

xtc_1001_silent: test
	./test -f ../ExampleData/M2A_S3_1001_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -m 10

xtc_1001_ring: test
	./test -f ../ExampleData/M2A_S3_1001_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -m 10 -r

run_silent: test
	./test -f ../ExampleData/M1A_S2_16.pdb -o ../ExampleData/test_out.pdb

supercell: test
	./test -f ../ExampleData/m1_tetra_102.xtc -t ../ExampleData/run1.top -o ../ExampleData/supercell.pdb -a

gmon.out: xtc_1001_silent

multithread:
	export OMP_NUM_THREADS=4

perf: gmon.out test
	gprof test gmon.out > perf_analysis.txt

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
	rm -f gmon.out perf_analysis.txt
	rm -f test
	rm -f ../ExampleData/test_out.pdb
	rm -f cluster_max_size.dat
	rm -f cluster_statistics.dat
	rm -f ring_statistics.dat