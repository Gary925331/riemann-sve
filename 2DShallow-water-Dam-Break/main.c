#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define NX 200          /* number of X cells */
#define NY 200          /* number of Y cells */
#define N (NX+2)*(NY+2)
#define NIF_X (NX+1)      /* number of X interfaces */
#define NIF_Y (NY+1)      /* number of Y interfaces */
#define NIF (NIF_X+2)*(NIF_Y+2)
#define ALPHA 1.0    /* advection speed */
#define L 1.0          /* domain length */
#define D 1.0          /* domain length */
#define DX (L / NX)    /* cell size */
#define DY (D / NY)    /* cell size */
#define MAX_TIMESTEPS 50000
#define T_FINAL 0.04
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
	for (int i = 0;i < NX+2;i++){
		for (int j = 0;j < NY+2;j++){
			int index = i*(NY+2)+j;
			if(i < (NX+2)/2){
				h[index] = 1;
			}else{
				h[index] = 1;
			}
//		printf("h[%d]=%f\n",i,h[i]);
			u[index] = 0;
			v[index] = 0;
			mass[index] = h[index];
			momentum_X[index] = h[index]*u[index];
			momentum_Y[index] = h[index]*v[index];

			mass_slope_X[index] = 0;
            		momentum_slope_X_X[index] = 0;
            		momentum_slope_X_Y[index] = 0;
            		mass_slope_Y[index] = 0;
            		momentum_slope_Y_X[index] = 0;
            		momentum_slope_Y_Y[index] = 0;
		//mass_slope[i] = 0;
		//momentum_slope[i] = 0;
		}
	}
	float time = 0;

	for (int timestep = 0; timestep < MAX_TIMESTEPS; timestep++){
		float Smax = 0.0;
		
		for (int i = 1;i < NX+1;i++){
			for (int j = 1;j < NY+1;j++){
				int index = i*(NY+2)+j;
				float S_L = fabs(u[index-NY-2]) + sqrt(g*h[index-NY-2]);
				float S_R = fabs(u[index]) + sqrt(g*h[index]);
				float S_B = fabs(v[index-1]) + sqrt(g*mass[index-1]);
                        	float S_T = fabs(v[index]) + sqrt(g*mass[index]);
				float S_local_max = S_L;
				if (S_R > S_local_max){
        				S_local_max = S_R;
    				}
    				if (S_B > S_local_max){
        				S_local_max = S_B;
    				}
    				if (S_T > S_local_max){
        				S_local_max = S_T;
    				}
    				if (S_local_max > Smax){
        				Smax = S_local_max;
    				}
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
		for(int i = 1;i < NX+1;i++){
			for(int j = 0;j < NY+2;j++){
				int index = i*(NY+2)+j;
				float forward = (mass[index+NY+2] - mass[index])/DX;
				float backward = (mass[index] - mass[index-NY-2])/DX;
				if(forward*backward < 0){
					mass_slope_X[index] = 0;
				}else{
					if(fabs(forward)<fabs(backward)){
						mass_slope_X[index] = forward;
					}else{
						mass_slope_X[index] = backward;
					}
				}
//			printf("mass_slope[%d] = %f\n",i,mass_slope[i]);
			}
		}
		for(int i = 0;i < NX+2;i++){
                        for(int j = 1;j < NY+1;j++){
                                int index = i*(NY+2)+j;
                                float forward = (mass[index+1] - mass[index])/DY;
                                float backward = (mass[index] - mass[index-1])/DY;
                                if(forward*backward < 0){
                                        mass_slope_Y[index] = 0;
                                }else{
                                        if(fabs(forward)<fabs(backward)){
                                                mass_slope_Y[index] = forward;
                                        }else{
                                                mass_slope_Y[index] = backward;
                                        }
                                }
//                      printf("mass_slope[%d] = %f\n",i,mass_slope[i]);
                        }
                }

//		printf("mass_slope[%d] = %f\n",i,mass_slope[i]);
		for(int i = 1;i < NX+1;i++){
			for(int j = 0;j < NY+2;j++){
				int index = i*(NY+2)+j;
                        	float forward = (momentum_X[index+NY+2] - momentum_X[index])/DX;
                        	float backward = (momentum_X[index] - momentum_X[index-NY-2])/DX;
                        	if(forward*backward < 0){
                                	momentum_slope_X_X[index] = 0;
                        	}else{
                                	if(fabs(forward)<fabs(backward)){
                                        	momentum_slope_X_X[index] = forward;
                                	}else{
                                        	momentum_slope_X_X[index] = backward;
                                	}
				}
			}
			//printf("momentum_slope[%d] = %f\n",i,momentum_slope[i]);
                }
		for(int i = 1;i < NX+1;i++){
			for(int j = 0;j < NY+2;j++){
                                int index = i*(NY+2)+j;
                        	float forward = (momentum_Y[index+NY+2] - momentum_Y[index])/DY;
                        	float backward = (momentum_Y[index] - momentum_Y[index-NY-2])/DY;
                        	if(forward*backward < 0){
                                	momentum_slope_X_Y[index] = 0;
                        	}else{
                                	if(fabs(forward)<fabs(backward)){
                                        	momentum_slope_X_Y[index] = forward;
                                	}else{
                                        	momentum_slope_X_Y[index] = backward;
                                	}
                        	}
			}
                        //printf("momentum_slope[%d] = %f\n",i,momentum_slope[i]);
                }
		for(int i = 0;i < NX+2;i++){
                        for(int j = 1;j < NY+1;j++){
                                int index = i*(NY+2)+j;
                                float forward = (momentum_X[index+1] - momentum_X[index])/DX;
                                float backward = (momentum_X[index] - momentum_X[index-1])/DX;
                                if(forward*backward < 0){
                                        momentum_slope_Y_X[index] = 0;
                                }else{
                                        if(fabs(forward)<fabs(backward)){
                                                momentum_slope_Y_X[index] = forward;
                                        }else{
                                                momentum_slope_Y_X[index] = backward;
                                        }
                        	}
			}
                        //printf("momentum_slope[%d] = %f\n",i,momentum_slope[i]);
                }
		for(int i = 1;i < NX+1;i++){
                        for(int j = 0;j < NY+2;j++){
                                int index = i*(NY+2)+j;
                                float forward = (momentum_Y[index+1] - momentum_Y[index])/DY;
                                float backward = (momentum_Y[index] - momentum_Y[index-1])/DY;
                                if(forward*backward < 0){
                                        momentum_slope_Y_Y[index] = 0;
                                }else{
                                        if(fabs(forward)<fabs(backward)){
                                                momentum_slope_Y_Y[index] = forward;
                                        }else{
                                                momentum_slope_Y_Y[index] = backward;
                                        }
                        	}
			}
                        //printf("momentum_slope[%d] = %f\n",i,momentum_slope[i]);
                }
		int wall = (NX+2)/2;
		for (int i = 1; i < NIF_X+2; i++){
			for (int j = 0; j < NY+2; j++){
				int index = i*(NY+2)+j;
				if(i == wall){
					float mass_l = mass[index-NY-2]; //+ 0.5*DX*mass_slope_X[index-NY-2];
					float momentum_X_l = momentum_X[index-NY-2]; //+ 0.5*DX*momentum_slope_X_X[index-NY-2];
					float momentum_Y_l = momentum_Y[index-NY-2]; //+ 0.5*DX*momentum_slope_X_Y[index-NY-2];
					float u_l = momentum_X_l / mass_l;
					float v_l = momentum_Y_l / mass_l;
					float mass_left = 0;
					float mom_left_X = momentum_X_l * u_l + 0.5*g*mass_l*mass_l;
					float mom_Left_Y = 0;

                                	mass_F[index] = mass_left;
                                	momentum_F_X[index] = mom_left_X;
                                	momentum_F_Y[index] = mom_Left_Y;
				}else{
            				float mass_l = mass[index-NY-2] + 0.5*DX*mass_slope_X[index-NY-2];
            				float mass_r = mass[index] - 0.5*DX*mass_slope_X[index];
            				float momentum_X_l = momentum_X[index-NY-2] + 0.5*DX*momentum_slope_X_X[index-NY-2];
            				float momentum_X_r = momentum_X[index] - 0.5*DX*momentum_slope_X_X[index];
					float momentum_Y_l = momentum_Y[index-NY-2] + 0.5*DX*momentum_slope_X_Y[index-NY-2];
                                	float momentum_Y_r = momentum_Y[index] - 0.5*DX*momentum_slope_X_Y[index];
            				float u_l = momentum_X_l / mass_l;
            				float u_r = momentum_X_r / mass_r;
					float v_l = momentum_Y_l / mass_l;
					float v_r = momentum_Y_r / mass_r;
            				float mass_left = mass_l * u_l;
            				float mass_right = mass_r * u_r;
            				float mom_left_X = momentum_X_l * u_l + 0.5*g*mass_l*mass_l;
            				float mom_right_X = momentum_X_r * u_r + 0.5*g*mass_r*mass_r;
					float mom_left_Y = momentum_Y_l * u_l;
                                	float mom_right_Y = momentum_Y_r * u_r;

            				float S_L = fabs(u_l) + sqrt(g*mass_l);
            				float S_R = fabs(u_r) + sqrt(g*mass_r);
					float S = S_L;
                                	if (S_R > S){
                                        	S = S_R;
                                	}

            				mass_F[index] = 0.5*(mass_left + mass_right) - 0.5*S*(mass_r - mass_l);
            				momentum_F_X[index] = 0.5*(mom_left_X + mom_right_X) - 0.5*S*(momentum_X_r - momentum_X_l);
					momentum_F_Y[index] = 0.5*(mom_left_Y + mom_right_Y) - 0.5*S*(momentum_Y_r - momentum_Y_l);
				}
			}
        	}
		for (int i=0; i<NY+2; i++){
			mass_F[i] = 0;
                	mass_F[N+i] = 0;
                	momentum_F_X[i] = 0.5*g*h[i]*h[i];
                	momentum_F_X[N+i] = 0.5*g*h[(NX+1)*(NY+2)+i]*h[(NX+1)*(NY+2)+i];
			momentum_F_Y[i] = 0;
                        momentum_F_Y[N+i] = 0;

		}
		for (int i = 0; i < NX+2; i++){
                        for (int j = 1; j < NIF_Y+2; j++){
                                int index = i*(NY+2)+j;
				//int index1 = i*(NIF_Y+2)+j;
                                float mass_B = mass[index-1] + 0.5*DY*mass_slope_Y[index-1];
                                float mass_T = mass[index] - 0.5*DY*mass_slope_Y[index];
                                float momentum_X_B = momentum_X[index-1] + 0.5*DY*momentum_slope_Y_X[index-1];
                                float momentum_X_T = momentum_X[index] - 0.5*DY*momentum_slope_Y_X[index];
				float momentum_Y_B = momentum_Y[index-1] + 0.5*DY*momentum_slope_Y_Y[index-1];
                                float momentum_Y_T = momentum_Y[index] - 0.5*DY*momentum_slope_Y_Y[index];
                                float u_B = momentum_X_B / mass_B;
                                float u_T = momentum_X_T / mass_T;
				float v_B = momentum_Y_B / mass_B;
				float v_T = momentum_Y_T / mass_T;
                                float mass_Bottom = mass_B * v_B;
                                float mass_Top = mass_T * v_T;
                                float mom_Bottom_X = momentum_X_B * v_B;
                                float mom_Top_X = momentum_X_T * v_T;
				float mom_Bottom_Y = momentum_Y_B * v_B + 0.5*g*mass_B*mass_B;
                                float mom_Top_Y = momentum_Y_T * v_T + 0.5*g*mass_T*mass_T;

                                float S_B = fabs(v_B) + sqrt(g*mass_B);
                                float S_T = fabs(v_T) + sqrt(g*mass_T);
                                float S = S_B;
                                if (S_T > S){
                                        S = S_T;
                                }

                                mass_G[index] = 0.5*(mass_Bottom + mass_Top) - 0.5*S*(mass_T - mass_B);
                                momentum_G_X[index] = 0.5*(mom_Bottom_X + mom_Top_X) - 0.5*S*(momentum_X_T - momentum_X_B);
				momentum_G_Y[index] = 0.5*(mom_Bottom_Y + mom_Top_Y) - 0.5*S*(momentum_Y_T - momentum_Y_B);
                        }
                }
		for (int i=0; i<NX+2; i++){
                        mass_G[i] = 0;
                        mass_G[N+i] = 0;
                        momentum_G_Y[i] = 0.5*g*h[i]*h[i];
                        momentum_G_Y[N+i] = 0.5*g*h[(NY+1)*(NX+2)+i]*h[(NY+1)*(NX+2)+i];
			momentum_G_X[i] = 0;
                        momentum_G_X[N+i] = 0;

                }
		// ghost mesh X direction
		for(int j = 0; j < NY+2; j++){
    			int left_ghost = 0*(NY+2) + j;
    			int left_inner = 1*(NY+2) + j;
    			int right_ghost = (NX+1)*(NY+2) + j;
    			int right_inner = NX*(NY+2) + j;
    
    			mass[left_ghost] = mass[left_inner];
    			momentum_X[left_ghost] = momentum_X[left_inner];
    			momentum_Y[left_ghost] = momentum_Y[left_inner];
    
    			mass[right_ghost] = mass[right_inner];
    			momentum_X[right_ghost] = momentum_X[right_inner];
    			momentum_Y[right_ghost] = momentum_Y[right_inner];
		}

		// ghost mesh Y direction
		for(int i = 0; i < NX+2; i++){
    			int bottom_ghost = i*(NY+2) + 0;
    			int bottom_inner = i*(NY+2) + 1;
    			int top_ghost = i*(NY+2) + (NY+1);
    			int top_inner = i*(NY+2) + NY;

    			mass[bottom_ghost] = mass[bottom_inner];
    			momentum_X[bottom_ghost] = momentum_X[bottom_inner];
    			momentum_Y[bottom_ghost] = momentum_Y[bottom_inner];
    
    			mass[top_ghost] = mass[top_inner];
    			momentum_X[top_ghost] = momentum_X[top_inner];
    			momentum_Y[top_ghost] = momentum_Y[top_inner];
		}
		
		for(int i = 1;i < NX+1;i++){
			for (int j = 1; j < NY+1; j++){
				int index = i*(NY+2)+j;
				//int index1 = i*(NIF_Y+2)+j;
				mass[index] = mass[index] - (DT*(mass_F[index+NY+2]-mass_F[index])/DX) - (DT*(mass_G[index+1]-mass_G[index])/DY);
				momentum_X[index] = momentum_X[index] - (DT*(momentum_F_X[index+NY+2]-momentum_F_X[index])/DX)-(DT*(momentum_G_X[index+1]-momentum_G_X[index])/DY);
				momentum_Y[index] = momentum_Y[index] - (DT*(momentum_F_Y[index+NY+2]-momentum_F_Y[index])/DX)-(DT*(momentum_G_Y[index+1]-momentum_G_Y[index])/DY);
//			printf("%f,%f\n",mass[i],momentum[i]);
			} 
		}
		for(int i = 1;i < NX+1;i++){
                        for (int j = 1; j < NY+1; j++){
                                int index = i*(NY+2)+j;
				h[index] = mass[index];
				u[index] = momentum_X[index]/mass[index];
				v[index] = momentum_Y[index]/mass[index];
			}
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
	float *momentum_slope_X_X;
	float *momentum_slope_X_Y;
	float *mass_slope_Y;
	float *momentum_slope_Y_X;
	float *momentum_slope_Y_Y;

        Allocate_memory(&u,&v,&mass_F,&momentum_F_X,&momentum_F_Y,&mass_G,&momentum_G_X,&momentum_G_Y,&mass,&momentum_X,&momentum_Y,&h,
	&mass_slope_X,&momentum_slope_X_X,&momentum_slope_X_Y,&mass_slope_Y,&momentum_slope_Y_X,&momentum_slope_Y_Y);
	Calculation(u,v,mass_F,momentum_F_X,momentum_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentum_Y,h,
	mass_slope_X,momentum_slope_X_X,momentum_slope_X_Y,mass_slope_Y,momentum_slope_Y_X,momentum_slope_Y_Y);
	FILE *fp = fopen("results.dat", "w");
    	for (int j = 1; j < NX+1; j++) {
		 for (int k = 1; k < NY+1; k++) {
			int index = j*(NY+2)+k;
			float X = (j+0.5)*DX;
			float Y = (k+0.5)*DY;
        		fprintf(fp, "%g\t%g\t%g\t%g\t%g\n", X,Y, mass[index],u[index],v[index]);
		}
    	}
    	fclose(fp);
        Free_memory(u,v,mass_F,momentum_F_X,momentum_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentum_Y,h,mass_slope_X,momentum_slope_X_X,momentum_slope_X_Y,mass_slope_Y,momentum_slope_Y_X,momentum_slope_Y_Y);
}

