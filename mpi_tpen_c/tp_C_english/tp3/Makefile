############################ -*- Mode: Makefile -*- ###########################
## Makefile --- MPI : TP3 : collectives communications
##
## Author          : Dimitri Lecas (CNRS/IDRIS) <dimitri.lecas@idris.fr>
## 
###############################################################################

# Compiles, options for compilation and link
#include ../arch/make_inc

OBJS  = pi.o

# Implicit rule for compilation
.SUFFIXES: .o .c
.c.o:
	$(CC) -c $(CFLAGS_TP3) $<

all: pi

pi: $(OBJS)
	$(CC) -o $@ $(LDFLAGS_TP3) $(OBJS)

.PHONY : clean exe

exe : pi
	mpiexec -n 6 ./pi
#SBATCH	sbatch pi.slurm

clean:
	rm -f $(OBJS) pi core
