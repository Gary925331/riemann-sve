#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define N 100          /* number of cells */
#define NIF (N+1)      /* number of interfaces */
#define ALPHA 1.0    /* advection speed */
#define L 1.0          /* domain length */
#define DX (L / N)    /* cell size */
#define MAX_TIMESTEPS 5000
#define T_FINAL 0.1
#define g 9.81
#define CFL 0.5

void Allocate_memory(float **u,float **mass_F,float **momentum_F,float **mass,float **momentum,float **h){
	*u = (float*)malloc(N*sizeof(float));
	*mass_F = (float*)malloc(NIF*sizeof(float));
	*momentum_F = (float*)malloc(NIF*sizeof(float));
	*mass = (float*)malloc(N*sizeof(float));
	*momentum = (float*)malloc(N*sizeof(float));
	*h = (float*)malloc(N*sizeof(float));
}
void Free_memory(float *u,float *mass_F,float *momentum_F,float *mass,float *momentum,float *h){
	free(u);
	free(mass_F);
	free(momentum_F);
	free(mass);
	free(momentum);
	free(h);
}
void Calcaulation(float *u,float *mass_F,float *momentum_F,float *mass,float *momentum,float *h){
	for (int i = 0;i < N;i++){
		if(i < N/2){
			h[i] = 10;
		}else{
			h[i] = 1;
		}
		printf("h[%d]=%f\n",i,h[i]);
		u[i] = 0;
		mass[i] = h[i];
		momentum[i] = h[i]*u[i];
	}
	float time = 0;
        float Smax;
	for (int timestep = 0; timestep < MAX_TIMESTEPS; timestep++){
		for (int i = 1;i < NIF-1;i++){
			float S_L = u[i-1] + sqrt(g*h[i-1]);
			float S_R = u[i] + sqrt(g*h[i]);
			float S_local_max;
			if (S_L > S_R){
				S_L = S_local_max;
			}else{
				S_R = S_local_max;
			}
			if(S_local_max > Smax){
				Smax = S_local_max;
			}
		}
		float DT = CFL*DX/Smax;
		time = time + DT;
        	if (time > T_FINAL) {
            		printf("Arrived at target time; stopping.\n");
            		break;
  		} else {
	            	printf("Ran out of timesteps before reaching target time.\n");
        	}
		for (int i = 1;i < NIF-1;i++){
			float mass_left = h[i-1]*u[i-1];
			float mass_right = h[i]*u[i];
			float momentum_left = h[i-1]*u[i-1]*u[i-1] + 0.5*g*h[i-1]*h[i-1];
			float momentum_right = h[i]*u[i]*u[i] + 0.5*g*h[i]*h[i];
                        mass_F[i] = 0.5*(mass_left + mass_right) + 0.5*Smax*(mass[i] - mass[i-1]);
			momentum_F[i] = 0.5*(momentum_left + momentum_right) +0.5*Smax*(momentum[i] - momentum[i-1]);
			mass_F[0] = mass_F[1];
			mass_F[NIF-1] = mass_F[NIF-2];
			momentum_F[0] = momentum_F[1];
			momentum_F[NIF-1] = momentum[NIF-2];
		}
		for(int i = 1;i < N;i++){
		}


	}
}
int main() {
        float *u;
        float *mass_F;
        float *momentum_F;
        float *mass;
        float *momentum;
        float *h;
        Allocate_memory(&u,&mass_F,&momentum_F,&mass,&momentum,&h);
	Calcaulation(u,mass_F,momentum_F,mass,momentum,h);
        Free_memory(u,mass_F,momentum_F,mass,momentum,h);
}

