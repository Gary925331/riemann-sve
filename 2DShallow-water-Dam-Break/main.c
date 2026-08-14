#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define NX 200          /* number of X cells */
#define NY 200          /* number of Y cells */
#define N NX*NY
#define NIF_X (NX+1)      /* number of X interfaces */
#define NIF_Y (NY+1)      /* number of Y interfaces */
#define NIF NIF_X*NIF_Y
#define ALPHA 1.0    /* advection speed */
#define L 1.0          /* domain length */
#define D 1.0          /* domain length */
#define DX (L / NX)    /* cell size */
#define DY (D / NY)    /* cell size */
#define MAX_TIMESTEPS 50000
#define T_FINAL 0.1
#define g 9.81
#define CFL 0.25

void Allocate_memory(float **u,float **v,float **mass_F,float **momentum_F_X,float **momentum_F_Y,float **mass_G,float **momentum_G_X,float **momentum_G_Y,float **mass,float **momentum_X,float **momentum_Y,float **h,float **mass_slope_X,float **momentum_slope_X_X,float **momentum_slope_X_Y,float **mass_slope_Y,float **momentum_slope_Y_X,float **momentum_slope_Y_Y){
	*u = (float*)malloc(N*sizeof(float));
	*v = (float*)malloc(N*sizeof(float));
	*mass_F = (float*)malloc(NIF*sizeof(float));
	*momentum_F_X = (float*)malloc(NIF*sizeof(float));
	*momentum_F_Y = (float*)malloc(NIF*sizeof(float));
	*mass_G = (float*)malloc(NIF*sizeof(float));
        *momentum_G_X = (float*)malloc(NIF*sizeof(float));
	*momentum_G_Y = (float*)malloc(NIF*sizeof(float));
	*mass = (float*)malloc(N*sizeof(float));
	*momentum_X = (float*)malloc(N*sizeof(float));
	*momentum_Y = (float*)malloc(N*sizeof(float));
	*h = (float*)malloc(N*sizeof(float));
	*mass_slope_X = (float*)malloc(N*sizeof(float));
	*momentum_slope_X_X = (float*)malloc(N*sizeof(float));
	*momentum_slope_X_Y = (float*)malloc(N*sizeof(float));
	*mass_slope_Y = (float*)malloc(N*sizeof(float));
        *momentum_slope_Y_X = (float*)malloc(N*sizeof(float));
	*momentum_slope_Y_Y = (float*)malloc(N*sizeof(float));
}
void Free_memory(float *u,float *v,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,float *momentum_Y,float *h,float *mass_slope_X,float *momentum_slope_X_X,float *momentum_slope_X_Y,float *mass_slope_Y,float *momentum_slope_Y_X,float *momentum_slope_Y_Y){
	free(u);
	free(v);
	free(mass_F);
	free(momentum_F_X);
	free(momentum_F_Y);
	free(mass_G);
        free(momentum_G_X);
	free(momentum_G_Y);
	free(mass);
	free(momentum_X);
	free(momentum_Y);
	free(h);
	free(mass_slope_X);
	free(momentum_slope_X_X);
	free(momentum_slope_X_Y);
	free(mass_slope_Y);
        free(momentum_slope_Y_X);
	free(momentum_slope_Y_Y);
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

void Calculation(float *u,float *v,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,float *momentum_Y,float *h,float *mass_slope_X,float *momentum_slope_X_X,float *momentum_slope_X_Y,float *mass_slope_Y,float *momentum_slope_Y_X,float *momentum_slope_Y_Y){
	for (int i = 0;i < N+2;i++){
		if(i < N/2){
			h[i] = 1;
		}else{
			h[i] = 0.1;
		}
//		printf("h[%d]=%f\n",i,h[i]);
		u[i] = 0;
		mass[i] = h[i];
		momentum_X[i] = h[i]*u[i];
		momentum_Y[i] = h[i]*v[i];
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
		for(int i = 1;i < NX-1;i++){
			for(int j = 0;j < NY;j++){
				int index = i*NY+j;
				float forward = (mass[index+NY] - mass[index])/DX;
				float backward = (mass[index] - mass[index-NY])/DX;
				if(forward*backward < 0){
					mass_slope_X[i] = 0;
				}else{
					if(fabs(forward)<fabs(backward)){
						mass_slope_X[i] = forward;
					}else{
						mass_slope_X[i] = backward;
					}
				}
//			printf("mass_slope[%d] = %f\n",i,mass_slope[i]);
			}
		}
//		printf("mass_slope[%d] = %f\n",i,mass_slope[i]);
		for(int i = 1;i < N;i++){
                        float forward = (momentum_X[i+1] - momentum_X[i])/DX;
                        float backward = (momentum_X[i] - momentum_X[i-1])/DX;
                        if(forward*backward < 0){
                                momentum_slope_X_X[i] = 0;
                        }else{
                                if(fabs(forward)<fabs(backward)){
                                        momentum_slope_X_X[i] = forward;
                                }else{
                                        momentum_slope_X_X[i] = backward;
                                }
			}
			//printf("momentum_slope[%d] = %f\n",i,momentum_slope[i]);
                }

		for (int i = 1; i < NIF_X-1; i++){
			for (int j = 0; j < NY; j++){
				int index = i*NY+j;
            			float mass_l = mass[index-NY] + 0.5*DX*mass_slope_X[index-NY];
            			float mass_r = mass[index] - 0.5*DX*mass_slope_X[index];
            			float momentum_X_l = momentum_X[index-NY] + 0.5*DX*momentum_slope_X_X[index-NY];
            			float momentum_X_r = momentum_X[index] - 0.5*DX*momentum_slope_X_X[index];
				float momentum_Y_l = momentum_Y[index-NY] + 0.5*DX*momentum_slope_X_Y[index-NY];
                                float momentum_Y_r = momentum_Y[index] - 0.5*DX*momentum_slope_X_Y[index];
            			float u_l = momentum_X_l / mass_l;
            			float u_r = momentum_X_r / mass_r;
				float v_l = momentum_Y_l / mass_l;
				float v_r = momentum_Y_r / mass_r;
            			float mass_left = mass_l * u_l;
            			float mass_right = mass_r * u_r;
            			float mom_left_X = momentum_X_l * u_l + 0.5*g*mass_l*mass_l;
            			float mom_right_X = momentum_X_r * u_r + 0.5*g*mass_r*mass_r;
				float mom_left_Y = momentum_Y_l * v_l;
                                float mom_right_Y = momentum_Y_r * v_r;

            			float S_L = u[i-1] + sqrt(g*mass[i-1]);
            			float S_R = u[i] + sqrt(g*mass[i]);
            			float S;
            			if (S_L > S_R) {
                			S = S_L;
            			} else {
                			S = S_R;
            			}
            			mass_F[i] = 0.5*(mass_left + mass_right) - 0.5*S*(mass_r - mass_l);
            			momentum_F_X[i] = 0.5*(mom_left_X + mom_right_X) - 0.5*S*(momentum_X_r - momentum_X_l);
				momentum_F_Y[i] = 0.5*(mom_left_Y + mom_right_Y) - 0.5*S*(momentum_Y_r - momentum_Y_l);
			}
        	}
		mass_F[0] = 0;
                mass_F[NIF-1] = 0;
                momentum_F_X[0] = 0.5*g*h[0]*h[0];
                momentum_F_Y[NIF-1] = 0.5*g*h[NIF-2]*h[NIF-2];
		for (int i = 1; i < NIF_X-1; i++){
                        for (int j = 0; j < NY-1; j++){
                                int index = i*NX+j;
                                float mass_B = mass[index-1] + 0.5*DX*mass_slope_Y[i-1];
                                float mass_T = mass[index] - 0.5*DX*mass_slope_Y[index];
                                float momentum_X_B = momentum_X[index-1] + 0.5*DX*momentum_slope[index-1];
                                float momentum_X_T = momentum_X[index] - 0.5*DX*momentum_slope[index];
				float momentum_Y_B = momentum_Y[index-1] + 0.5*DX*momentum_slope[index-1];
                                float momentum_Y_T = momentum_Y[index] - 0.5*DX*momentum_slope[index];
                                float u_B = momentum_X_B / mass_B;
                                float u_T = momentum_X_T / mass_T;
				float v_B = momentum_Y_B / mass_B;
				float v_T = momentum_Y_T / mass_T;
                                float mass_Bottom = mass_B * u_B;
                                float mass_Top = mass_T * u_T;
                                float mom_Bottom_X = momentum_X_B * v_B;
                                float mom_Top_X = momentum_X_T * v_T;
				float mom_Bottom_Y = momentum_Y_B * u_B + 0.5*g*mass_B*mass_B;
                                float mom_Top_Y = momentum_Y_T * u_T + 0.5*g*mass_T*mass_T;

                                float S_L = u[index-1] + sqrt(g*mass[index-1]);
                                float S_R = u[index] + sqrt(g*mass[index]);
                                float S;
                                if (S_L > S_R) {
                                        S = S_L;
                                } else {
                                        S = S_R;
                                }
                                mass_G[index] = 0.5*(mass_Bottom + mass_Top) - 0.5*S*(mass_T - mass_B);
                                momentum_G_X[index] = 0.5*(mom_Bottom_X + mom_Top_X) - 0.5*S*(momentum_X_T - momentum_X_B);
				momentum_G_Y[index] = 0.5*(mom_Bottom_Y + mom_Top_Y) - 0.5*S*(momentum_Y_T - momentum_Y_B);
                        }
                }

		for(int i = 0;i < N;i++){
			mass[index] = mass[index] - DT*(mass_F[index+1]-mass_F[index])/DX;
			momentum[index] = momentum[index] - DT*(momentum_F[index+1]-momentum_F[index])/DX;
//			printf("%f,%f\n",mass[i],momentum[i]); 
		}
		for(int i = 0;i < N;i++){
			h[index] = mass[index];
			u[index] = momentum[index]/mass[index];
		}
		time = time + DT;

	}
}
int main() {
        float *u;
	float *v;
        float *mass_F;
        float *momentum_F_X;
	float *momentum_F_Y;
	float *mass_G;
        float *momentum_G_X;
        float *momentum_G_Y;
        float *mass;
        float *momentum_X;
	float *momentum_Y;
        float *h;
	float *mass_slope_X;
	float *momentum_slope_X;
	float *mass_slope_Y;
        float *momentum_slope_Y;
        Allocate_memory(&u,&v,&mass_F,&momentum_F_X,&momentum_F_Y,&mass_G,&momentum_G_X,&momentum_G_Y,&mass,&momentum_X,&momentum_Y,&h,&mass_slope_X,&momentum_slope_X,&mass_slope_Y,&momentum_slope_Y);
	Calculation(u,v,mass_F,momentum_F_X,momentim_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentim_Y,h,mass_slope_X,momentum_slope_X_X,momentum_slope_X_Y,mass_slope_Y,momentum_slope_Y_X,momentum_slope_Y_Y);
	FILE *fp = fopen("results.dat", "w");
    	for (int j = 0; j < N; j++) {
		float X = (j+0.5)*DX;
        	fprintf(fp, "%g\t%.15e\t%g\n", X, mass[j],u[j]);
    	}
    	fclose(fp);
        Free_memory(u,v,mass_F,momentum_F_X,momentim_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentim_Y,h,mass_slope_X,momentum_slope_X,mass_slope_Y,momentum_slope_Y);
}

