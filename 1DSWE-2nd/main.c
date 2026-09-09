#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define N 1000          /* number of cells */
#define NIF (N+1)      /* number of interfaces */
#define ALPHA 1.0    /* advection speed */
#define L 1.0          /* domain length */
#define DX (L / N)    /* cell size */
#define MAX_TIMESTEPS 50000
#define T_FINAL 0.1
#define g 9.81
#define CFL 0.25

void Allocate_memory(float **u,float **mass_F,float **momentum_F,float **mass,float **momentum,float **h,float **mass_slope,float **momentum_slope){
	*u = (float*)malloc(N*sizeof(float));
	*mass_F = (float*)malloc(NIF*sizeof(float));
	*momentum_F = (float*)malloc(NIF*sizeof(float));
	*mass = (float*)malloc(N*sizeof(float));
	*momentum = (float*)malloc(N*sizeof(float));
	*h = (float*)malloc(N*sizeof(float));
	*mass_slope = (float*)malloc(N*sizeof(float));
	*momentum_slope = (float*)malloc(N*sizeof(float));
}
void Free_memory(float *u,float *mass_F,float *momentum_F,float *mass,float *momentum,float *h,float *mass_slope,float *momentum_slope){
	free(u);
	free(mass_F);
	free(momentum_F);
	free(mass);
	free(momentum);
	free(h);
	free(mass_slope);
	free(momentum_slope);
}
/*float MINMOD(float QL_rho, float QC_rho, float QR_rho, float dx){
	float dU_dx;
	float Forward = (QR_rho - QC_rho) / dx;
	float Backward = (QC_rho - QL_rho) / dx;
	if (Backward * Forward < 0){
		dU_dx = 0;
	} else if ( fabs(Forward) < fabs(Backward) ){
			dU_dx = Forward;
		} else {
			dU_dx = Backward;
	}
	return dU_dx;
}*/
float Compare_S(float u1,float u2,float h1,float h2){
	float S_local_max;
        float S_L = u1 + sqrt(g*h1);
        float S_R = u2 + sqrt(g*h2);
        if (S_L > S_R){
        	S_local_max = S_L;
        }else{
                S_local_max = S_R;
        }
	return S_local_max;
}

void Calculation(float *u,float *mass_F,float *momentum_F,float *mass,float *momentum,float *h,float *mass_slope,float *momentum_slope){
	for (int i = 0;i < N+2;i++){
		if(i < N/2){
			h[i] = 1;
		}else{
			h[i] = 0.1;
		}
//		printf("h[%d]=%f\n",i,h[i]);
		u[i] = 0;
		mass[i] = h[i];
		momentum[i] = h[i]*u[i];
		//mass_slope[i] = 0;
		//momentum_slope[i] = 0;
	}
	float time = 0;

	for (int timestep = 0; timestep < MAX_TIMESTEPS; timestep++){
		float Smax = 0.0;
		float S_local_max = 0.0;
		for (int i = 1;i < NIF-1;i++){
			float S_L = u[i-1] + sqrt(g*h[i-1]);
			float S_R = u[i] + sqrt(g*h[i]);
			if (S_L > S_R){
				S_local_max = S_L;
			}else{
				S_local_max = S_R;
			}
			if(S_local_max > Smax){
				Smax = S_local_max;
			}
		}
		float DT = CFL*DX/Smax;
		printf("%f\n",Smax);
//		time = time + DT;
        	if (time > T_FINAL) {
            		//printf("Arrived at target time; stopping.\n");
            		break;
  		} else {
	            	//printf("Ran out of timesteps before reaching target time.\n");
   		}
		for(int i = 1;i < N-1;i++){
			float forward = (mass[i+1] - mass[i])/DX;
			float backward = (mass[i] - mass[i-1])/DX;
			if(forward*backward < 0){
				mass_slope[i] = 0;
			}else{
				if(fabs(forward)<fabs(backward)){
					mass_slope[i] = forward;
				}else{
					mass_slope[i] = backward;
				}
			}
//			printf("mass_slope[%d] = %f\n",i,mass_slope[i]);
		}
		mass_slope[0]=0;
                mass_slope[N-1]=0;
//		printf("mass_slope[%d] = %f\n",i,mass_slope[i]);
		for(int i = 1;i < N;i++){
                        float forward = (momentum[i+1] - momentum[i])/DX;
                        float backward = (momentum[i] - momentum[i-1])/DX;
                        if(forward*backward < 0){
                                momentum_slope[i] = 0;
                        }else{
                                if(fabs(forward)<fabs(backward)){
                                        momentum_slope[i] = forward;
                                }else{
                                        momentum_slope[i] = backward;
                                }
			}
			//printf("momentum_slope[%d] = %f\n",i,momentum_slope[i]);
                }
		momentum_slope[0]=0;
		momentum_slope[N-1]=0;

		/*for (int i = 1;i < NIF-1;i++){
			float mass_l = mass[i-1] + 0.5*DX*mass_slope[i-1];
                        float mass_r = mass[i] - 0.5*DX*mass_slope[i];
                        float momentum_l = momentum[i-1] + 0.5*DX*momentum_slope[i-1];
                        float momentum_r = momentum[i] - 0.5*DX*momentum_slope[i];
			float mass_left = mass[i-1]*u[i-1];
			float mass_right = mass[i]*u[i];
			float momentum_left = h[i-1]*u[i-1]*u[i-1] + 0.5*g*h[i-1]*h[i-1];
			float momentum_right = h[i]*u[i]*u[i] + 0.5*g*h[i]*h[i];
			float S = Compare_S(u[i-1],u[i],h[i-1],h[i]);

                        mass_F[i] = 0.5*(mass_left + mass_right) - 0.5*S*(mass_r - mass_l);
			momentum_F[i] = 0.5*(momentum_left + momentum_right) - 0.5*S*(momentum_r - momentum_l);
		}*/
		for (int i = 1; i < NIF-1; i++){
            		float mass_l = mass[i-1] + 0.5*DX*mass_slope[i-1];
            		float mass_r = mass[i] - 0.5*DX*mass_slope[i];
            		float momentum_l = momentum[i-1] + 0.5*DX*momentum_slope[i-1];
            		float momentum_r = momentum[i] - 0.5*DX*momentum_slope[i];
            		float u_l = momentum_l / mass_l;
            		float u_r = momentum_r / mass_r;
            		float mass_left = mass_l * u_l;
            		float mass_right = mass_r * u_r;
            		float mom_left = momentum_l * u_l + 0.5*g*mass_l*mass_l;
            		float mom_right = momentum_r * u_r + 0.5*g*mass_r*mass_r;

            		float S_L = u[i-1] + sqrt(g*mass[i-1]);
            		float S_R = u[i] + sqrt(g*mass[i]);
            		float S;
            		if (S_L > S_R) {
                		S = S_L;
            		} else {
                		S = S_R;
            		}
            		mass_F[i] = 0.5*(mass_left + mass_right) - 0.5*S*(mass_r - mass_l);
            		momentum_F[i] = 0.5*(mom_left + mom_right) - 0.5*S*(momentum_r - momentum_l);
        	}
		mass_F[0] = 0;
                mass_F[NIF-1] = 0;
                momentum_F[0] = 0.5*g*h[0]*h[0];
                momentum_F[NIF-1] = 0.5*g*h[NIF-2]*h[NIF-2];

		for(int i = 0;i < N;i++){
			mass[i] = mass[i] - DT*(mass_F[i+1]-mass_F[i])/DX;
			momentum[i] = momentum[i] - DT*(momentum_F[i+1]-momentum_F[i])/DX;
//			printf("%f,%f\n",mass[i],momentum[i]); 
		}
		for(int i = 0;i < N;i++){
			h[i] = mass[i];
			u[i] = momentum[i]/mass[i];
		}
		time = time + DT;

	}
}
int main() {
        float *u;
        float *mass_F;
        float *momentum_F;
        float *mass;
        float *momentum;
        float *h;
	float *mass_slope;
	float *momentum_slope;
        Allocate_memory(&u,&mass_F,&momentum_F,&mass,&momentum,&h,&mass_slope,&momentum_slope);
	Calculation(u,mass_F,momentum_F,mass,momentum,h,mass_slope,momentum_slope);
	FILE *fp = fopen("results.dat", "w");
    	for (int j = 0; j < N; j++) {
		float X = (j+0.5)*DX;
        	fprintf(fp, "%g\t%.15e\t%g\n", X, mass[j],u[j]);
    	}
    	fclose(fp);
        Free_memory(u,mass_F,momentum_F,mass,momentum,h,mass_slope,momentum_slope);
}

