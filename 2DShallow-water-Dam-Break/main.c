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
#define L 200.0          /* domain length */
#define D 200.0          /* domain length */
#define DX (L / NX)    /* cell size */
#define DY (D / NY)    /* cell size */
#define MAX_TIMESTEPS 50000
#define T_FINAL 7.2
#define g 9.81
#define CFL 0.5

void Allocate_memory(float **u,float **v,float **s,float **mass_F,float **momentum_F_X,float **momentum_F_Y,float **mass_G,float **momentum_G_X,float **momentum_G_Y,float **mass,float **momentum_X,float **momentum_Y,float **h,float **mass_slope_X,float **momentum_slope_X_X,float **momentum_slope_X_Y,float **mass_slope_Y,float **momentum_slope_Y_X,float **momentum_slope_Y_Y){
	*u = (float*)malloc(N*sizeof(float));
	*v = (float*)malloc(N*sizeof(float));
	*s = (float*)malloc(N*sizeof(float));
	*mass_F = (float*)malloc(NIF*sizeof(float));
	*momentum_F_X = (float*)malloc(NIF*sizeof(float));
	*momentum_F_Y = (float*)malloc(NIF*sizeof(float));
	*mass_G = (float*)malloc(NIF*sizeof(float));
        *momentum_G_X = (float*)malloc(NIF*sizeof(float));
	*momentum_G_Y = (float*)malloc(NIF*sizeof(float));
	*mass = (float*)malloc(NIF*sizeof(float));
	*momentum_X = (float*)malloc(NIF*sizeof(float));
	*momentum_Y = (float*)malloc(NIF*sizeof(float));
	*h = (float*)malloc(N*sizeof(float));
	*mass_slope_X = (float*)malloc(N*sizeof(float));
	*momentum_slope_X_X = (float*)malloc(N*sizeof(float));
	*momentum_slope_X_Y = (float*)malloc(N*sizeof(float));
	*mass_slope_Y = (float*)malloc(N*sizeof(float));
        *momentum_slope_Y_X = (float*)malloc(N*sizeof(float));
	*momentum_slope_Y_Y = (float*)malloc(N*sizeof(float));
}
void Free_memory(float *u,float *v,float *s,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,float *momentum_Y,float *h,float *mass_slope_X,float *momentum_slope_X_X,float *momentum_slope_X_Y,float *mass_slope_Y,float *momentum_slope_Y_X,float *momentum_slope_Y_Y){
	free(u);
	free(v);
	free(s);
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

void Calculation(float *u,float *v,float *s,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,float *momentum_Y,float *h,float *mass_slope_X,float *momentum_slope_X_X,float *momentum_slope_X_Y,float *mass_slope_Y,float *momentum_slope_Y_X,float *momentum_slope_Y_Y){
	for (int i = 0;i < NX+2;i++){
		for (int j = 0;j < NY+2;j++){
			int index = i*(NY+2)+j;
			if(i < (NX+2)/2){
				h[index] = 10;
			}else{
				h[index] = 5;
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
		}
	}
	float time = 0;

	for (int timestep = 0; timestep < MAX_TIMESTEPS; timestep++){
		float Smax_X = 0.0;
		float Smax_Y = 0.0;
		for (int i = 1;i < NX+1;i++){
			for (int j = 1;j < NY+1;j++){
				int index = i*(NY+2)+j;
				float S_R = fabs(u[index]) + sqrt(g*h[index]);
				if (S_R > Smax_X) {
            				Smax_X = S_R;
        			}
                        	float S_T = fabs(v[index]) + sqrt(g*mass[index]);
				if (S_T > Smax_Y) {
            				Smax_Y = S_T;
        			}
			}
		}
		float term_X = Smax_X / DX;
		float term_Y = Smax_Y / DY;
		float max_term;
		if (term_X > term_Y) {
    			max_term = term_X;
		} else {
    			max_term = term_Y;
		}

		float DT = CFL / max_term;
		printf("%f\n",max_term);
//		time = time + DT;
        	if (time > T_FINAL) {
            		//printf("Arrived at target time; stopping.\n");
            		break;
  		} else {
	            	//printf("Ran out of timesteps before reaching target time.\n");
   		}

		for(int j = 0; j < NY+2; j++){
            		mass[0*(NY+2)+j] = mass[1*(NY+2)+j];
            		momentum_X[0*(NY+2)+j] = momentum_X[1*(NY+2)+j]; 
            		momentum_Y[0*(NY+2)+j] = momentum_Y[1*(NY+2)+j];
            
            		mass[(NX+1)*(NY+2)+j] = mass[NX*(NY+2)+j];
            		momentum_X[(NX+1)*(NY+2)+j] = momentum_X[NX*(NY+2)+j]; 
            		momentum_Y[(NX+1)*(NY+2)+j] = momentum_Y[NX*(NY+2)+j];

		}
        	for(int i = 0; i < NX+2; i++){
            		mass[i*(NY+2)+0] = mass[i*(NY+2)+1];
            		momentum_X[i*(NY+2)+0] = momentum_X[i*(NY+2)+1];
            		momentum_Y[i*(NY+2)+0] = -momentum_Y[i*(NY+2)+1]; // 下牆反彈
            
            		mass[i*(NY+2)+NY+1] = mass[i*(NY+2)+NY];
            		momentum_X[i*(NY+2)+NY+1] = momentum_X[i*(NY+2)+NY];
            		momentum_Y[i*(NY+2)+NY+1] = -momentum_Y[i*(NY+2)+NY]; // 上牆反彈
        	}
		
		for(int i = 0; i < NX+2; i++){
                        for (int j = 0; j < NY+2; j++){
                                int index = i*(NY+2)+j;
                                h[index] = mass[index];
                                u[index] = momentum_X[index]/mass[index];
                                v[index] = momentum_Y[index]/mass[index];
                        }
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
                        	float forward = (u[index+NY+2] - u[index])/DX;
                        	float backward = (u[index] - u[index-NY-2])/DX;
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
                        	float forward = (v[index+NY+2] - v[index])/DX;
                        	float backward = (v[index] - v[index-NY-2])/DX;
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
                                float forward = (u[index+1] - u[index])/DY;
                                float backward = (u[index] - u[index-1])/DY;
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
                                float forward = (v[index+1] - v[index])/DY;
                                float backward = (v[index] - v[index-1])/DY;
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
		//X direction flux
		for (int i = 1; i < NIF_X+2; i++){
			for (int j = 0; j < NY+2; j++){
				int index = i*(NY+2)+j;
            				float mass_l = mass[index-NY-2] ;//+ 0.5*DX*mass_slope_X[index-NY-2];
            				float mass_r = mass[index] ;//- 0.5*DX*mass_slope_X[index];
            				float u_X_l = u[index-NY-2] ;//+ 0.5*DX*momentum_slope_X_X[index-NY-2];
            				float u_X_r = u[index] ;//- 0.5*DX*momentum_slope_X_X[index];
					float v_Y_l = v[index-NY-2] ;//+ 0.5*DX*momentum_slope_X_Y[index-NY-2];
                                	float v_Y_r = v[index] ;//- 0.5*DX*momentum_slope_X_Y[index];
					if (j <= 96 || j >= 171) {
            					if (i == 101) {
                					mass_r = mass_l;
                					u_X_r = -u_X_l;
                					v_Y_r = v_Y_l;
            					}
            					else if (i == 102) {
                					mass_l = mass_r;
                					u_X_l = -u_X_r;
                					v_Y_l = v_Y_r;
            					}
        				}
            				float u_l = u_X_l;
            				float u_r = u_X_r;
					float v_l = v_Y_l;
					float v_r = v_Y_r;
            				float mass_left = mass_l * u_l;
            				float mass_right = mass_r * u_r;
            				float mom_left_X = mass_l * u_l * u_l + 0.5*g*mass_l*mass_l;
            				float mom_right_X = mass_r * u_r * u_r + 0.5*g*mass_r*mass_r;
					float mom_left_Y = mass_l * v_l * u_l;
                                	float mom_right_Y = mass_r * v_r * u_r; 

            				float S_L = fabs(u_l) + sqrt(g*mass_l);
            				float S_R = fabs(u_r) + sqrt(g*mass_r);
                                	float S;
					if (S_L > S_R) {
    						S = S_L;
					} else {
    						S = S_R;
					}
            				mass_F[index] = 0.5*(mass_left + mass_right) - 0.5*S*(mass_r - mass_l);
            				momentum_F_X[index] = 0.5*(mom_left_X + mom_right_X) - 0.5*S*(mass_r * u_r - mass_l * u_l);
					momentum_F_Y[index] = 0.5*(mom_left_Y + mom_right_Y) - 0.5*S*(mass_r * v_r - mass_l * v_l);
			}
        	}

		//Y direction flux
		for (int i = 0; i < NX+2; i++){
                        for (int j = 1; j < NIF_Y+2; j++){
                                int index = i*(NY+2)+j;
				//int index1 = i*(NIF_Y+2)+j;
                                float mass_B = mass[index-1] ;//+ 0.5*DY*mass_slope_Y[index-1];
                                float mass_T = mass[index] ;//- 0.5*DY*mass_slope_Y[index];
                                float u_X_B = u[index-1] ;//+ 0.5*DY*momentum_slope_Y_X[index-1];
                                float u_X_T = u[index] ;//- 0.5*DY*momentum_slope_Y_X[index];
				float v_Y_B = v[index-1] ;//+ 0.5*DY*momentum_slope_Y_Y[index-1];
                                float v_Y_T = v[index] ;//- 0.5*DY*momentum_slope_Y_Y[index];
				if(i == 101){
					if(j == 97 ){
						mass_B = mass_T;
						u_X_B = u_X_T;
						v_Y_B = -v_Y_T;
					}
					if(j == 171){
						mass_T = mass_B;
						u_X_T = u_X_B;
						v_Y_T = -v_Y_B;
					}
				}
                                float u_B = u_X_B;
                                float u_T = u_X_T;
				float v_B = v_Y_B;
				float v_T = v_Y_T;
                                float mass_Bottom = mass_B * v_B;
                                float mass_Top = mass_T * v_T;
                                float mom_Bottom_X = mass_B * u_B * v_B;
                                float mom_Top_X = mass_T * u_T * v_T;
				float mom_Bottom_Y = mass_B * v_B * v_B + 0.5*g*mass_B*mass_B;
                                float mom_Top_Y = mass_T * v_T * v_T + 0.5*g*mass_T*mass_T;

                                float S_B = fabs(v_B) + sqrt(g*mass_B);
                                float S_T = fabs(v_T) + sqrt(g*mass_T);
                                float S;
                                if (S_T > S_B){
                                        S = S_T;
                                }else{
					S = S_B;
				}

                                mass_G[index] = 0.5*(mass_Bottom + mass_Top) - 0.5*S*(mass_T - mass_B);
                                momentum_G_X[index] = 0.5*(mom_Bottom_X + mom_Top_X) - 0.5*S*(mass_T * u_T - mass_B * u_B);
				momentum_G_Y[index] = 0.5*(mom_Bottom_Y + mom_Top_Y) - 0.5*S*(mass_T * v_T - mass_B * v_B);
                        }
                }
		for(int i = 1;i < NX+1;i++){
			for (int j = 1; j < NY+1; j++){
				int index = i*(NY+2)+j;
				//int index1 = i*(NIF_Y+2)+j;
				if (i == 101 && (j <= 96 || j >= 171)) {
            				continue; 
        			}
				mass[index] = mass[index] - (DT*(mass_F[index+NY+2]-mass_F[index])/DX) - (DT*(mass_G[index+1]-mass_G[index])/DY);
				momentum_X[index] = momentum_X[index] - (DT*(momentum_F_X[index+NY+2]-momentum_F_X[index])/DX)-(DT*(momentum_G_X[index+1]-momentum_G_X[index])/DY);
				momentum_Y[index] = momentum_Y[index] - (DT*(momentum_F_Y[index+NY+2]-momentum_F_Y[index])/DX)-(DT*(momentum_G_Y[index+1]-momentum_G_Y[index])/DY);
//			printf("%f,%f\n",mass[i],momentum[i]);
			} 
		}
		for(int i = 0;i < NX+2;i++){
                        for (int j = 0; j < NY+2; j++){
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
	float *s;
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

        Allocate_memory(&u,&v,&s,&mass_F,&momentum_F_X,&momentum_F_Y,&mass_G,&momentum_G_X,&momentum_G_Y,&mass,&momentum_X,&momentum_Y,&h,
	&mass_slope_X,&momentum_slope_X_X,&momentum_slope_X_Y,&mass_slope_Y,&momentum_slope_Y_X,&momentum_slope_Y_Y);
	Calculation(u,v,s,mass_F,momentum_F_X,momentum_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentum_Y,h,
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
        Free_memory(u,v,s,mass_F,momentum_F_X,momentum_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentum_Y,h,mass_slope_X,momentum_slope_X_X,momentum_slope_X_Y,mass_slope_Y,momentum_slope_Y_X,momentum_slope_Y_Y);
}

