all: main main_musl

main: test3.c
	gcc -O3 -march=native -ffast-math test3.c -c -lm
	gcc -O3 -march=native -ffast-math test3.c -o main -lm

main_musl: test3.c
	musl-gcc -O3 -march=native -ffast-math test3.c -o main_musl -lm

clean:
	rm -f main main_musl
