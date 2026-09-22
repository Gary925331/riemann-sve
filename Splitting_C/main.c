#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "memory.h"
#include "flux.h"
#include "initial.h"
#include "minmod.h"
#include "state.h"
#include "ghost.h"

#define NX 1000          /* number of X cells */
#define NY 1000          /* number of Y cells */
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
#define CFL 0.02

void time_calculation(float *u,float *v,float *s,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,
float *momentum_Y,float *h,float *h_slope_X,float *u_slope_X_X,float *v_slope_X_Y,float *h_slope_Y,float *u_slope_Y_X,float *v_slope_Y_Y,int type){
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
		ghost(u,v,s,mass_F,momentum_F_X,momentum_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentum_Y,h,
        	h_slope_X,u_slope_X_X,v_slope_X_Y,h_slope_Y,u_slope_Y_X,v_slope_Y_Y,NX,NY);

        	minmod(u,v,s,mass_F,momentum_F_X,momentum_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentum_Y,h,
        	h_slope_X,u_slope_X_X,v_slope_X_Y,h_slope_Y,u_slope_Y_X,v_slope_Y_Y,NX,NY,DX,DY);

        	Calculation(u,v,s,mass_F,momentum_F_X,momentum_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentum_Y,h,
        	h_slope_X,u_slope_X_X,v_slope_X_Y,h_slope_Y,u_slope_Y_X,v_slope_Y_Y,type,NX,NY,DX,DY,g);

        	state(u,v,s,mass_F,momentum_F_X,momentum_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentum_Y,h,
        	h_slope_X,u_slope_X_X,v_slope_X_Y,h_slope_Y,u_slope_Y_X,v_slope_Y_Y,NX,NY,DX,DY,DT);

		#pragma omp master
		{
		time = time + DT;
		}
		#pragma omp barrier

	}
	}//end of parallel
}
int main(int argc, char *argv[]) {
	int num_cores = 8; 
    	if (argc > 1) {
        	num_cores = atoi(argv[1]);
    	}
	int type = FLUX_HLL;
	if (argc > 2) {
    		type = atoi(argv[2]);   // ./program 8 0  → HLL；./program 8 1 → LF
	}
	if (type == FLUX_HLL) {
    		printf(">>> 目前使用的通量算法：HLL (type=%d)\n", type);
	} else if (type == FLUX_RUS) {
    		printf(">>> 目前使用的通量算法：Rusnaov (type=%d)\n", type);
	} else {
    		printf(">>> 警告：type=%d 不是有效值，程式會走 else 分支\n", type);
	}
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
	float *h_slope_X;
	float *u_slope_X_X;
	float *v_slope_X_Y;
	float *h_slope_Y;
	float *u_slope_Y_X;
	float *v_slope_Y_Y;
	
        Allocate_memory(&u,&v,&s,&mass_F,&momentum_F_X,&momentum_F_Y,&mass_G,&momentum_G_X,&momentum_G_Y,&mass,&momentum_X,&momentum_Y,&h,
	&h_slope_X,&u_slope_X_X,&v_slope_X_Y,&h_slope_Y,&u_slope_Y_X,&v_slope_Y_Y,N,NIF);
	
	time_t start_date;
    	time(&start_date);
    	printf("Simulation started at: %s", ctime(&start_date));

    	double start_wtime = omp_get_wtime(); // 取得開始的精確秒數
	
	Initial(u,v,s,mass_F,momentum_F_X,momentum_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentum_Y,h,
        h_slope_X,u_slope_X_X,v_slope_X_Y,h_slope_Y,u_slope_Y_X,v_slope_Y_Y,NX,NY);
	
	omp_set_num_threads(num_cores);

	time_calculation(u,v,s,mass_F,momentum_F_X,momentum_F_Y,mass_G,momentum_G_X,momentum_G_Y,mass,momentum_X,momentum_Y,h,
        h_slope_X,u_slope_X_X,v_slope_X_Y,h_slope_Y,u_slope_Y_X,v_slope_Y_Y,type);

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
        Free_memory(&u,&v,&s,&mass_F,&momentum_F_X,&momentum_F_Y,&mass_G,&momentum_G_X,&momentum_G_Y,&mass,&momentum_X,&momentum_Y,&h,
        &h_slope_X,&u_slope_X_X,&v_slope_X_Y,&h_slope_Y,&u_slope_Y_X,&v_slope_Y_Y);

}
