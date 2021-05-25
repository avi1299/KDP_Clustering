IDIR =./include
SDIR =./src
CC=gcc
CFLAGS=-I$(IDIR) -pg

ODIR=./obj
LDIR =./lib

LIBS=-lm

_DEPS = DFS_clustering.h HPO.h input.h stack.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = HPO.o stack.o DFS_clustering.o  input.o  test.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

test: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

run: test
	./test -f M1C_S2.pdb

run_silent: test
	./test -f M1C_S2.pdb > /dev/null

gmon.out: run_silent

perf: gmon.out test
	gprof test gmon.out > perf_analysis.txt

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
	rm -f gmon.out perf_analysis.txt
	rm -f test