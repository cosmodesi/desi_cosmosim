CC = gcc 
MPICC = mpicc
CFLAGS = -O3 -Wall -std=c99 -D_POSIX_C_SOURCE=200809L
LIBS = -fopenmp
INCL =
OPTS = $(INCL) $(LIBS)
SRCS_OMP = jobfork.c cmdlist.c subprocess.c
SRCS_MPI = $(SRCS_OMP) mpiengine.c
EXEC_OMP = jobfork_omp
EXEC_MPI = jobfork_mpi

all: omp mpi

omp:
	$(CC) $(CFLAGS) -o $(EXEC_OMP) $(SRCS_OMP) $(OPTS) -DCMD_OMP

mpi:
	$(MPICC) $(CFLAGS) -o $(EXEC_MPI) $(SRCS_MPI) $(OPTS) -DCMD_MPI

clean:
	rm -f $(EXEC_OMP) $(EXEC_MPI)
