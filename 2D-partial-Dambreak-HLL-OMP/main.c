#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
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
	float Smax_X;
        float Smax_Y;
	float DT;
	int stop_flag = 0;
	#pragma omp parallel
	{
	for (int timestep = 0; timestep < MAX_TIMESTEPS; timestep++){
		#pragma omp single
		{
			Smax_X = 0.0;
			Smax_Y = 0.0;
		}
		#pragma omp for reduction(max:Smax_X, Smax_Y)
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
		#pragma omp master
                {
		float term_X = Smax_X / DX;
		float term_Y = Smax_Y / DY;
		float max_term;
		if (term_X > term_Y) {
    			max_term = term_X;
		} else {
    			max_term = term_Y;
		}

		DT = CFL / max_term;
		
//		printf("%f\n",max_term);
//		time = time + DT;
        	if (time > T_FINAL) {
            		//printf("Arrived at target time; stopping.\n");
            		stop_flag = 1;
  		} 
		}
		#pragma omp barrier
		if (stop_flag == 1) {
            		break;
        	}
		#pragma omp for
		for(int j = 0; j < NY+2; j++){
            		mass[0*(NY+2)+j] = mass[1*(NY+2)+j];
            		momentum_X[0*(NY+2)+j] = -momentum_X[1*(NY+2)+j]; // 左牆反彈
            		momentum_Y[0*(NY+2)+j] = momentum_Y[1*(NY+2)+j];
            
            		mass[(NX+1)*(NY+2)+j] = mass[NX*(NY+2)+j];
            		momentum_X[(NX+1)*(NY+2)+j] = -momentum_X[NX*(NY+2)+j]; // 右牆反彈
            		momentum_Y[(NX+1)*(NY+2)+j] = momentum_Y[NX*(NY+2)+j];

		}
		#pragma omp for
        	for(int i = 0; i < NX+2; i++){
            		mass[i*(NY+2)+0] = mass[i*(NY+2)+1];
            		momentum_X[i*(NY+2)+0] = momentum_X[i*(NY+2)+1];
            		momentum_Y[i*(NY+2)+0] = -momentum_Y[i*(NY+2)+1]; // 下牆反彈
            
            		mass[i*(NY+2)+NY+1] = mass[i*(NY+2)+NY];
            		momentum_X[i*(NY+2)+NY+1] = momentum_X[i*(NY+2)+NY];
            		momentum_Y[i*(NY+2)+NY+1] = -momentum_Y[i*(NY+2)+NY]; // 上牆反彈
        	}
		#pragma omp for
		for(int i = 0; i < NX+2; i++){
                        for (int j = 0; j < NY+2; j++){
                                int index = i*(NY+2)+j;
                                h[index] = mass[index];
                                u[index] = momentum_X[index]/mass[index];
                                v[index] = momentum_Y[index]/mass[index];
                        }
                }
		#pragma omp for
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
		#pragma omp for
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
		#pragma omp for
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
		#pragma omp for
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
		#pragma omp for
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
		#pragma omp for
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
		#pragma omp for
		//X direction flux
		for (int i = 1; i < NX+2; i++){
			for (int j = 0; j < NY+2; j++){
				int index = i*(NY+2)+j;
            				float mass_l = mass[index-NY-2] + 0.5*DX*mass_slope_X[index-NY-2];
            				float mass_r = mass[index] - 0.5*DX*mass_slope_X[index];
            				float u_X_l = u[index-NY-2] + 0.5*DX*momentum_slope_X_X[index-NY-2];
            				float u_X_r = u[index] - 0.5*DX*momentum_slope_X_X[index];
					float v_Y_l = v[index-NY-2] + 0.5*DX*momentum_slope_X_Y[index-NY-2];
                                	float v_Y_r = v[index] - 0.5*DX*momentum_slope_X_Y[index];
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

            				float SL1 = (u_l) - sqrt(g*mass_l);
            				float SR1 = (u_r) - sqrt(g*mass_r);
                                	float S1; //SL
					if (SL1 > SR1) {
    						S1 = SR1;
					} else {
    						S1 = SL1;
					}
					float SL2 = (u_l) + sqrt(g*mass_l);
					float SR2 = (u_r) + sqrt(g*mass_r);
					float S2; //SR
					if (SL2 > SR2) {
                                                S2 = SL2;
                                        } else {
                                                S2 = SR2;
                                        }
					if(S1 > 0){
						mass_F[index] = mass_left;
						momentum_F_X[index] = mom_left_X;
						momentum_F_Y[index] = mom_left_Y;
					}else if(S2 < 0){
						mass_F[index] = mass_right;
                                                momentum_F_X[index] = mom_right_X;
                                                momentum_F_Y[index] = mom_right_Y;
					}else{
            					mass_F[index] = (S2*mass_left - S1*mass_right)/(S2-S1) + S1*S2*(mass_r - mass_l)/(S2-S1);
            					momentum_F_X[index] = (S2*mom_left_X - S1*mom_right_X)/(S2-S1) + S1*S2*(mass_r * u_r - mass_l * u_l)/(S2-S1);
						momentum_F_Y[index] = (S2*mom_left_Y - S1*mom_right_Y)/(S2-S1) + S1*S2*(mass_r * v_r - mass_l * v_l)/(S2-S1);
					}
			}
        	}
		#pragma omp for
		//Y direction flux
		for (int i = 0; i < NX+2; i++){
                        for (int j = 1; j < NY+2; j++){
                                int index = i*(NY+2)+j;
				//int index1 = i*(NIF_Y+2)+j;
                                float mass_B = mass[index-1] + 0.5*DY*mass_slope_Y[index-1];
                                float mass_T = mass[index] - 0.5*DY*mass_slope_Y[index];
                                float u_X_B = u[index-1] + 0.5*DY*momentum_slope_Y_X[index-1];
                                float u_X_T = u[index] - 0.5*DY*momentum_slope_Y_X[index];
				float v_Y_B = v[index-1] + 0.5*DY*momentum_slope_Y_Y[index-1];
                                float v_Y_T = v[index] - 0.5*DY*momentum_slope_Y_Y[index];
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

                                float SB1 = (v_B) - sqrt(g*mass_B);
                                float ST1 = (v_T) - sqrt(g*mass_T);
                                float S1; //SL
                                if (ST1 > SB1){
                                        S1 = SB1;
                                }else{
					S1 = ST1;
				}
				float SB2 = (v_B) + sqrt(g*mass_B);
                                float ST2 = (v_T) + sqrt(g*mass_T);
                                float S2; //SR
                                if (ST2 > SB2){
                                        S2 = ST2;
                                }else{
                                        S2 = SB2;
                                }
				if(S1 > 0){
					mass_G[index] = mass_Bottom;
					momentum_G_X[index] = mom_Bottom_X;
					momentum_G_Y[index] = mom_Bottom_Y;
				}else if(S2 < 0){
					mass_G[index] = mass_Top;
                                        momentum_G_X[index] = mom_Top_X;
                                        momentum_G_Y[index] = mom_Top_Y;
				}else{
                                	mass_G[index] = (S2*mass_Bottom - S1*mass_Top)/(S2-S1) + S1*S2*(mass_T - mass_B)/(S2-S1);
                                	momentum_G_X[index] = (S2*mom_Bottom_X - S1*mom_Top_X)/(S2-S1) + S1*S2*(mass_T * u_T - mass_B * u_B)/(S2-S1);
					momentum_G_Y[index] = (S2*mom_Bottom_Y - S1*mom_Top_Y)/(S2-S1) + S1*S2*(mass_T * v_T - mass_B * v_B)/(S2-S1);
				}
                        }
                }
		#pragma omp for
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
		#pragma omp for
		for(int i = 0;i < NX+2;i++){
                        for (int j = 0; j < NY+2; j++){
                                int index = i*(NY+2)+j;
				h[index] = mass[index];
				u[index] = momentum_X[index]/mass[index];
				v[index] = momentum_Y[index]/mass[index];
			}
		}
		#pragma omp master
		{
		time = time + DT;
		}
		#pragma omp barrier
	}
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
	time_t start_date;
    	time(&start_date);
    	printf("Simulation started at: %s", ctime(&start_date));

    	double start_wtime = omp_get_wtime(); // 取得開始的精確秒數
	omp_set_num_threads(8);

	Calculation(u,v,s,mass_F,momentum_F_X,momentum_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentum_Y,h,
	mass_slope_X,momentum_slope_X_X,momentum_slope_X_Y,mass_slope_Y,momentum_slope_Y_X,momentum_slope_Y_Y);
	
	double end_wtime = omp_get_wtime(); // 取得結束的精確秒數

    	time_t end_date;
    	time(&end_date);
    	printf("Simulation finished at: %s", ctime(&end_date));

    	// 印出精確的執行時間
    	printf("Total Execution Time: %f seconds\n", end_wtime - start_wtime);

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

