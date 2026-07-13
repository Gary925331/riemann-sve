#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define N 100          /* number of cells */
#define NIF (N+1)      /* number of interfaces */
#define ALPHA 1.0    /* advection speed */
#define L 1.0          /* domain length */
#define DX (L / N)    /* cell size */

void Allocate_memory(float *u,float *mass_F,float *momentum_F,float *mass,float *momentum,float *h){
	u = (float*)malloc(N*sizeof(float));
	mass_F = (float*)malloc(NIF*sizeof(float));
	momentum_F = (float*)malloc(NIF*sizeof(float));
	mass = (float*)malloc(N*sizeof(float));
	momentum = (float*)malloc(N*sizeof(float));
	h = (float*)malloc(N*sizeof(float));
}
void Free_memory(float *u,float *mass_F,float *momentum_F,float *mass,float *momentum,float *h){
	free(u);
	free(mass_F);
	free(momentum_F);
	free(mass);
	free(momentum);
	free(h);
}
int main() {
        float *u;
        float *mass_F;
        float *momentum_F;
        float *mass;
        float *momentum;
        float *h;
        Allocate_memory(u,mass_F,momentum_F,mass,momentum,h);
        Free_memory(u,mass_F,momentum_F,mass,momentum,h);
        return 0;
}

