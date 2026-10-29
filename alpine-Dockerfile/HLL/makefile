SRC    = main.c flux.c ghost.c initial.c memory.c minmod.c state.c mclimiter.c 
HDR    = $(wildcard *.h)
CFLAGS = -O3 -ffast-math

all: main main_musl

main: $(SRC) $(HDR)
	gcc $(CLFAGS) $(SRC) -c main.c -lm
	gcc $(CFLAGS) $(SRC) -o main -lm

main_musl: $(SRC) $(HDR)
	musl-gcc $(CFLAGS) $(SRC) -o main_musl -lm

clean:
	rm -f main main_musl results_*cores_run*.dat results.dat

.PHONY: all clean
