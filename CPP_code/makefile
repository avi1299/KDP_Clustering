IDIR =./include
SDIR =./src
CC=g++
CFLAGS=-I$(IDIR) -pg -DGMX_DOUBLE=0 -fopenmp

ODIR=./obj
LDIR =./lib

LIBS=-lgromacs  -lm

_DEPS = DFS_clustering.h HPO.h input.h stack.h output.h graph.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = HPO.o stack.o DFS_clustering.o  input.o  output.o graph.o test.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

test: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -L$(LDIR) $(LIBS)

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