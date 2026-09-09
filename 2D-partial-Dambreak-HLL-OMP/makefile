all:
	gcc -fopenmp -O3 main.c -o main -lm

test: all
	@for cores in 1 2 4 8 16; do \
		echo "========================================"; \
		echo ">> Testing with $$cores cores..."; \
		./main $$cores; \
		mv results.dat results_$${cores}_cores.dat; \
	done
	@echo "All tests completed!"
