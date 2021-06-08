IDIR =./include
SDIR =./src
ODIR=./obj
LDIR =./lib
CC=g++
CFLAGS=-I$(IDIR) -pg -DGMX_DOUBLE=0 -fopenmp
MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
MKFILE_DIR := $(dir $(MKFILE_PATH))
LD_LIBRARY_PATH=$(MKFILE_DIR)lib

LIBS=-lgromacs  -lm

_DEPS = DFS_clustering.h HPO.h input.h stack.h output.h graph.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = HPO.o stack.o DFS_clustering.o  input.o  output.o graph.o test.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))



$(ODIR)/%.o:$(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -L$(LDIR) $(LIBS) -Wl,-rpath=$(LD_LIBRARY_PATH)

test:$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -L$(LDIR) $(LIBS) -Wl,-rpath=$(LD_LIBRARY_PATH)

run: test
	./test -f ../Datafiles/M1A_S2_16.pdb -o ../Datafiles/out.pdb -v1

run_silent: test
	./test -f ../Datafiles/M1A_S2_16.pdb -o ../Datafiles/out.pdb

gmon.out: run_silent


perf: gmon.out test
	gprof test gmon.out > perf_analysis.txt

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
	rm -f gmon.out perf_analysis.txt
	rm -f test