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

clusterIt:$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -L$(LDIR) $(LIBS) -Wl,-rpath=$(LD_LIBRARY_PATH)

run: clusterIt
	./clusterIt -f ../ExampleData/M1A_S2_16.pdb -o ../ExampleData/test_out.pdb -v1

small: clusterIt
	./clusterIt -f ../ExampleData/M1A_S1_Avi.xtc -t ../ExampleData/M1A_S1.top -o ../ExampleData/test_out.pdb -m 100

#xtc_1: clusterIt
#	./clusterIt -f ../ExampleData/M2A_S3_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -v2

#xtc_1_silent: clusterIt
#	./clusterIt -f ../ExampleData/M2A_S3_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -m 10

xtc_1001: clusterIt
	./clusterIt -f ../ExampleData/M2A_S3_1001_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -v2

xtc_1001_silent: clusterIt
	./clusterIt -f ../ExampleData/M2A_S3_1001_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -m 10

xtc_1001_ring: clusterIt
	./clusterIt -f ../ExampleData/M2A_S3_1001_Avi.xtc -t ../ExampleData/M2A_S3.top -o ../ExampleData/test_out.pdb -m 10 -r

run_silent: clusterIt
	./clusterIt -f ../ExampleData/M1A_S2_16.pdb -o ../ExampleData/test_out.pdb

supercell: clusterIt
	./clusterIt -f ../ExampleData/m1_tetra_102.xtc -t ../ExampleData/run1.top -o ../ExampleData/supercell.pdb -a

minicell_pdb: clusterIt
	./clusterIt -f ../ExampleData/avi_M1.pdb -o ../ExampleData/minicell_pdb.pdb

minicell_xtc: clusterIt
	./clusterIt -f ../ExampleData/avi_M1.xtc -t ../ExampleData/run2.top -o ../ExampleData/minicell_xtc.pdb -r


bsh: clusterIt
	./clusterIt -f ../ExampleData/M1B_Sh_avi.xtc -t ../ExampleData/run156.top -o ../ExampleData/M1B_Sh_out.pdb -m 10 

bsh_r: clusterIt
	./clusterIt -f ../ExampleData/M1B_Sh_avi.xtc -t ../ExampleData/run156.top -r -b -d

bs4_r: clusterIt
	./clusterIt -f ../ExampleData/M1B_S4_avi.xtc -t ../ExampleData/run1664.top -r -b -d

bs1: clusterIt
	./clusterIt -f ../ExampleData/M1B_S1_avi.xtc -t ../ExampleData/run312.top -o ../ExampleData/M1B_S1_out.pdb -m 10 

bs2: clusterIt
	./clusterIt -f ../ExampleData/M1B_S2_avi.xtc -t ../ExampleData/run624.top -o ../ExampleData/M1B_S2_out.pdb -m 10 

bs3: clusterIt
	./clusterIt -f ../ExampleData/M1B_S3_avi.xtc -t ../ExampleData/run1248.top -o ../ExampleData/M1B_S3_out.pdb -m 10 

bs4: clusterIt
	./clusterIt -f ../ExampleData/M1B_S4_avi.xtc -t ../ExampleData/run1664.top -o ../ExampleData/M1B_S4_out.pdb -m 10 

csh: clusterIt
	./clusterIt -f ../ExampleData/M1C_Sh_avi.xtc -t ../ExampleData/run156.top -o ../ExampleData/M1C_Sh_out.pdb -m

cs4: clusterIt
	./clusterIt -f ../ExampleData/M1C_S4_avi.xtc -t ../ExampleData/run1664.top -o ../ExampleData/M1C_S4_out.pdb -m 10


bs2pbc: clusterIt
	./clusterIt -f ../ExampleData/M1B_S2_pbc_avi.xtc -t ../ExampleData/run624.top -o ../ExampleData/M1B_S2_pbc_out.pdb -m 10





gmon.out: xtc_1001_silent

multithread:
	export OMP_NUM_THREADS=4

perf: gmon.out clusterIt
	gprof clusterIt gmon.out > perf_analysis.txt

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
	rm -f gmon.out perf_analysis.txt
	rm -f clusterIt
	rm -f ../ExampleData/test_out.pdb
	rm -f cluster_max_size.dat
	rm -f cluster_statistics.dat
	rm -f ring_statistics.dat
	rm -f Kstatistics.dat
