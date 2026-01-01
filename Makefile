all: runseq runomp runmpi

seq: main.c 
	gcc main.c -O2 -lm -o build/firefly

omp: openmp.c 
	gcc openmp.c -O2 -fopenmp -lm -o build/firefly_omp

mpi: mpi.c 
	mpicc mpi.c -O2 -lm -o build/firefly_mpi

runseq: seq
	/usr/bin/time -f "Czas: %e s" ./build/firefly

runomp: omp
	OMP_NUM_THREADS=8 /usr/bin/time -f "Czas: %e s" ./build/firefly_omp

runmpi: mpi
	/usr/bin/time -f "Czas: %e s" mpirun -np 4 ./build/firefly_mpi
