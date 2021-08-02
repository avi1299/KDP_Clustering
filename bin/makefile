IDIR =./include
SDIR =./src
ODIR=./obj
LDIR =./lib
CC=g++
CFLAGS=-I$(IDIR) -DGMX_DOUBLE=0 -fopenmp -mcmodel=medium #-pg -g
MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
MKFILE_DIR := $(dir $(MKFILE_PATH))
LD_LIBRARY_PATH=$(MKFILE_DIR)lib

LIBS=-lgromacs  -lm

_DEPS = DFS_clustering.h HPO.h input.h stack.h output.h graph.h constants.h charge.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = HPO.o stack.o DFS_clustering.o  input.o  output.o graph.o charge.o test.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))



$(ODIR)/%.o:$(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -L$(LDIR) $(LIBS) -Wl,-rpath=$(LD_LIBRARY_PATH)

test:$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -L$(LDIR) $(LIBS) -Wl,-rpath=$(LD_LIBRARY_PATH)

run: test
	./test -f ../ExampleData/M1A_S2_16.pdb -o ../ExampleData/test_out.pdb -v1

#xtc_1: test
#	./test -f ../ExampleData/M2A_S3_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -v2

#xtc_1_silent: test
#	./test -f ../ExampleData/M2A_S3_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -m 10

xtc_1001: test
	./test -f ../ExampleData/M2A_S3_1001_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -v2

xtc_1001_silent: test
	./test -f ../ExampleData/M2A_S3_1001_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -m 10

xtc_1001_silent_strong: test
	./test -f ../ExampleData/M2A_S3_1001_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -m 10 -r

run_silent: test
	./test -f ../ExampleData/M1A_S2_16.pdb -o ../ExampleData/test_out.pdb

test_out: test
	./test -f ../ExampleData/m1_tetra_102.xtc -t ../ExampleData/model1.top -o ../ExampleData/test_out102.pdb -m 10

gmon.out: xtc_1001_silent

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