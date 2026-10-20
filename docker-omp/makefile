all:
	gcc -fopenmp -O3 main.c flux.c ghost.c initial.c memory.c minmod.c state.c mclimiter.c -o main -lm

clean:
	rm -f main results_*cores_run*.dat results.dat

