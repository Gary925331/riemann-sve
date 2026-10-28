all: main main_musl

main: test.c
	gcc -fopenmp -O3 test.c -o main -lm

main_musl: test.c
	musl-gcc -O3 test.c -o main_musl -lm

clean:
	rm -f main main_musl
