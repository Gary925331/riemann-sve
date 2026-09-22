all:
	gcc -fopenmp -O3 main.c flux.c ghost.c initial.c memory.c minmod.c state.c -o main -lm

bench: all
	@for cores in 1 2 4; do \
		for run in 1 2 3; do \
			echo "=================================="; \
			echo "用 $$cores 核心跑，第 $$run 次"; \
			echo "=================================="; \
			./main $$cores 0; \
			mv results.dat results_$${cores}cores_run$${run}.dat; \
			echo "已存成 results_$${cores}cores_run$${run}.dat"; \
			echo ""; \
		done; \
	done
	@echo "全部跑完，產生的檔案："
	@ls -la results_*cores_run*.dat

clean:
	rm -f main results_*cores_run*.dat results.dat

