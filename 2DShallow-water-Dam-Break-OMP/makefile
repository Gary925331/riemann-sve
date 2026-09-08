all:
	gcc -fopenmp -O3 main.c -o main -lm

test: all
	@for cores in 1 2 4; do \
		echo "========================================"; \
		echo ">> Testing with $$cores cores..."; \
		for run in 1 2 3 4 5; do\
			echo "  - Run $$run:";\
			./main $$cores | grep "Total Execution Time";\
		done;\
	done
	@echo "All tests completed!"
