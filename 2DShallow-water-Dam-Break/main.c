
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
		// ghost mesh X direction
                for(int j = 1; j < NY+1; j++){
			int idx_R = 1*(NY+2)+j;      // 內部網格 Real mesh
			int idx_F_left = 1*(NY+2)+j; // interface
                        
			float h_R = mass[idx_R]+0.5*DX*mass_slope_X[idx_R];
			float momX_R = momentum_X[idx_R] - 0.5*DX*momentum_slope_X_X[idx_R];
            		float momY_R = momentum_Y[idx_R] - 0.5*DX*momentum_slope_X_Y[idx_R];
            		float u_R = momX_R / h_R;
            		float v_R = momY_R / h_R;
			//ghost 
			float h_G = h_R;
			float u_G = -u_R;
			float v_G = v_R;
			//U
			float F_R_mass = h_R * u_R;
            		float F_R_momX = h_R * u_R * u_R + 0.5 * g * h_R * h_R;
            		float F_R_momY = h_R * u_R * v_R;
            		float F_G_mass = h_G * u_G;
            		float F_G_momX = h_G * u_G * u_G + 0.5 * g * h_G * h_G;
            		float F_G_momY = h_G * u_G * v_G;
			float S = fabs(u_R) + sqrt(g * h_R); 
			
			mass_F[idx_F_left] = 0.5 * (F_G_mass + F_R_mass) - 0.5*S*(h_R - h_G);
			momentum_F_X[idx_F_left] = 0.5*(F_G_momX + F_R_momX) - 0.5*S*(momentum_X[idx_R] - (h_G*u_G));
			momentum_F_Y[idx_F_left] = 0.5*(F_G_momY + F_R_momY) - 0.5*S*(momentum_Y[idx_R] - (h_G*v_G));
			
			int idx_L = NX*(NY+2)+j;         // Real mesh
			int idx_F_right = (NX+1)*(NY+2)+j; // interface i=NX+1
            		float h_L = mass[idx_L]+0.5*DX*mass_slope_X[idx_L];
			float momX_L = momentum_X[idx_L] + 0.5*DX*momentum_slope_X_X[idx_L];
            		float momY_L = momentum_Y[idx_L] + 0.5*DX*momentum_slope_X_Y[idx_L];
            		float u_L = momX_L / h_L;
            		float v_L = momY_L / h_L;
			float h_G_r = h_L;
            		float u_G_r = -u_L;
            		float v_G_r = v_L;
            		float F_L_mass = h_L * u_L;
            		float F_L_momX = h_L * u_L * u_L + 0.5 * g * h_L * h_L;
            		float F_L_momY = h_L * u_L * v_L;
            		float F_G_r_mass = h_G_r * u_G_r;
            		float F_G_r_momX = h_G_r * u_G_r * u_G_r + 0.5 * g * h_G_r * h_G_r;
            		float F_G_r_momY = h_G_r * u_G_r * v_G_r;
            		float S1 = fabs(u_L) + sqrt(g * h_L);
            		mass_F[idx_F_right] = 0.5 * (F_L_mass + F_G_r_mass) - 0.5 * S1 * (h_G_r - h_L);
            		momentum_F_X[idx_F_right] = 0.5 * (F_L_momX + F_G_r_momX) - 0.5 * S1 * ((h_G_r*u_G_r) - momentum_X[idx_L]);
            		momentum_F_Y[idx_F_right] = 0.5 * (F_L_momY + F_G_r_momY) - 0.5 * S1 * ((h_G_r*v_G_r) - momentum_Y[idx_L]);   
                }

                // ghost mesh Y direction
                for(int i = 1; i < NX+1; i++){
                        int idx_T = i*(NY+2)+1;      // Real mesh
            		int idx_G_bot = i*(NY+2)+1;  //interface j=1
            		float h_T = mass[idx_T]+0.5*DY*mass_slope_Y[idx_T];
            		float momX_T = momentum_X[idx_T] - 0.5*DY*momentum_slope_Y_X[idx_T];
            		float momY_T = momentum_Y[idx_T] - 0.5*DY*momentum_slope_Y_Y[idx_T];
            		float u_T = momX_T / h_T;
            		float v_T = momY_T / h_T;
            		//Ghost
            		float h_G = h_T;
            		float u_G = u_T;
            		float v_G = -v_T;
            		float G_T_mass = h_T * v_T;
            		float G_T_momX = h_T * u_T * v_T;
            		float G_T_momY = h_T * v_T * v_T + 0.5 * g * h_T * h_T;
            		float G_G_mass = h_G * v_G;
            		float G_G_momX = h_G * u_G * v_G;
            		float G_G_momY = h_G * v_G * v_G + 0.5 * g * h_G * h_G;
            		float S_bot = fabs(v_T) + sqrt(g * h_T);
            		mass_G[idx_G_bot] = 0.5 * (G_G_mass + G_T_mass) - 0.5 * S_bot * (h_T - h_G);
            		momentum_G_X[idx_G_bot] = 0.5 * (G_G_momX + G_T_momX) - 0.5 * S_bot * (momentum_X[idx_T] - (h_G*u_G));
            		momentum_G_Y[idx_G_bot] = 0.5 * (G_G_momY + G_T_momY) - 0.5 * S_bot * (momentum_Y[idx_T] - (h_G*v_G));
         
            		int idx_B = i*(NY+2)+NY;       // Real mesh
            		int idx_G_top = i*(NY+2)+NY+1; // interface j=NY+1
            		float h_B = mass[idx_B]+0.5*DY*mass_slope_Y[idx_B];
            		float momX_B = momentum_X[idx_B] + 0.5*DY*momentum_slope_Y_X[idx_B];
            		float momY_B = momentum_Y[idx_B] + 0.5*DY*momentum_slope_Y_Y[idx_B];
            		float u_B = momX_B / h_B;
            		float v_B = momY_B / h_B;
            		//Ghost
            		float h_G_t = h_B;
            		float u_G_t = u_B;
            		float v_G_t = -v_B;
            		float G_B_mass = h_B * v_B;
            		float G_B_momX = h_B * u_B * v_B;
            		float G_B_momY = h_B * v_B * v_B + 0.5 * g * h_B * h_B;
            		float G_G_t_mass = h_G_t * v_G_t;
            		float G_G_t_momX = h_G_t * u_G_t * v_G_t;
            		float G_G_t_momY = h_G_t * v_G_t * v_G_t + 0.5 * g * h_G_t * h_G_t;
            		float S_top = fabs(v_B) + sqrt(g * h_B);
            		mass_G[idx_G_top] = 0.5 * (G_B_mass + G_G_t_mass) - 0.5 * S_top * (h_G_t - h_B);
            		momentum_G_X[idx_G_top] = 0.5 * (G_B_momX + G_G_t_momX) - 0.5 * S_top * ((h_G_t*u_G_t) - momentum_X[idx_B]);
            		momentum_G_Y[idx_G_top] = 0.5 * (G_B_momY + G_G_t_momY) - 0.5 * S_top * ((h_G_t*v_G_t) - momentum_Y[idx_B]);
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
		//X direction flux
		for (int i = 1; i < NIF_X+2; i++){
			for (int j = 0; j < NY+2; j++){
				int index = i*(NY+2)+j;
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

		//Y direction flux
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

