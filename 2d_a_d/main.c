#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NX 800          /* number of X cells */
#define NY 800
#define N (NX*NY)
#define NIF_X (NX+1)      /* number of interfaces */
#define NIF_Y (NY+1)
#define U 0.5    /* advection speed X */
#define V 0.25   /* advection speed Y */
#define L 1.0          /* domain length */
#define H 1.0
#define DX (L/NX)    /* cell size */
#define DY (H/NY)
#define DT 0.0001 /* time step size */   
#define T_FINAL 1     /* final time */
#define MAX_TIMESTEPS 10000
#define NP 8
#define alpha 0.000024 //diffusion speed
#define DEBUG 1

void Compute_Fluxes(const float *T, float *F,float *W,float *D,float time)
{
    /* compute flux at all interfaces j = 0..n (n+1 interfaces)
       interface j is between left = (j-1) and right = j (mod n) */

    // Only compute the inside interfaces; leave the outside 2 interfaces for later.
//	#pragma omp parallel
//	{
	#pragma omp for   
    	for (int j = 1; j < NIF_X-1; j++) {
		for (int k = 0; k < NY; k++){
			int index1 = j*NY+k;
			float Left_F = U*T[index1-NY];
			float Right_F = U*T[index1];
			float right_T = T[index1];
			float left_T = T[index1-NY]; 

			F[index1] = 0.5 * (Left_F + Right_F) - 0.25*(right_T - left_T);
        	}
	}
	#pragma omp for
	for (int i=0; i<NY; i++){
                F[i]=F[i+NY];
                F[NX*NY+i]=F[(NX-1)*NY+i];
	}
	#pragma omp for
	for (int j = 0; j < NX; j++) {
                for (int k = 1; k < NIF_Y; k++){
                        int index1 = j*NIF_Y+k;
			int index = j*NY+k;
                        float Top_W = V*T[index];
                        float Bottom_W = V*T[index-1];
                        float Top_T = T[index]; 
                        float Bottom_T = T[index-1]; 

                        W[index1] = 0.5 * (Top_W + Bottom_W) - 0.25*(Top_T - Bottom_T);
                }
        }
//	printf("step %g Interface %g\n",time/DT,F[4050]);
	#pragma omp for
        for (int i=0; i<NX; i++){
                W[i*NIF_Y]=W[i*NIF_Y+1];
                W[(i+1)*NIF_Y-1]=W[(i+1)*NIF_Y-2];
	}
	#pragma omp for
	for (int j = 0; j < NX; j++) {
                for (int k = 0; k < NY; k++){
                        int index = j*NY+k;
                        float Bottom,Top,Left,Right;
                        if(k == 0){
                                Bottom = T[index];
                        }else{
                                Bottom = T[index-1];
                        }
                        if(k == (NY-1)){
                                Top = T[index];
                        }else{
                                Top = T[index+1];
                        }
			if(j == 0){
                                Left = T[index];
                        }else{
                                Left = T[index-NY];
                        }
                        if(j == (NX-1)){
                                Right = T[index];
                        }else{
                                Right = T[index+NY];
                        }

                        D[index] = (alpha/DY/DY)*(Left+Right+Top+Bottom-4*T[index]);
		}
	}
/*	for (int j = 1; j < NIF_Y; j++) {
                for (int k = 0; k < NX; k++){
                        int index1 = j*NX+k;
                        W[index1] = (-alpha)*((T[index1] - T[index1-NX]) / DY);
                }
        }
	printf("step %g Interface %g\n",time/DT,W[3060]);
	for (int i=0; i<NX; i++){
                W[i]=W[i+NX];
                W[NX*NY+i]=W[NX*(NY-1)+i];
        }*/


       		//central flux
       		// F[j] = 0.5 * (left_F + right_F);


       		// Rusanov flux
       		// F[j] = 0.5 * (left_F + right_F) - 0.25*(right_u - left_u);


//	}//end of parallel
    	// Set dF/dx = 0 on left and right ends of our domain
}


void Update_State(const float *F, float *T,float *W,float *D,float *Tnew) {

    /*
    Solving
        dT/dt + dF/dx = 0
    So
        T* = T - dt*dF/dx
    */
//	#pragma omp parallel
//	{
	#pragma omp for
    	for (int j = 0; j < NX; j++) {
                for (int k = 0; k < NY; k++){
                        int index = j*NY+k;
			int index1 = j*NIF_Y+k;
        		Tnew[index] = T[index] - ((DT/DX)*(F[index+NY] - F[index])) - ((DT/DY)*(W[index1+1]-W[index1])) + DT*D[index];
    		}
	}
	#pragma omp for
	for (int j = 0; j < NX; j++) {
                for (int k = 0; k < NY; k++){
                        int index1 = j*NY+k; 
                        T[index1] = Tnew[index1];
                }
        }

}



int main(void)
{
	float *T;
	float *F;
	float *A;
	float *Tnew;
	float *W;
	float *D;
	T = (float*)malloc(NIF_X*NIF_Y*sizeof(float));
	Tnew = (float*)malloc(NIF_X*NIF_Y*sizeof(float));
	F = (float*)malloc(NIF_X*NY*sizeof(float));
	A = (float*)malloc(NIF_X*NIF_Y*sizeof(float));
	W = (float*)malloc(NIF_X*NIF_Y*sizeof(float));
	D = (float*)malloc(NIF_X*NIF_Y*sizeof(float));

    	if (!T || !F) {
        	fprintf(stderr, "allocation failed\n");
        	return 1;
    	}

    /* initial condition */
	omp_set_num_threads(NP);
	float Total_error = 0.0;

	FILE *pFile;
        if (DEBUG) printf("Saving results\n");
        pFile = fopen("results.txt", "w");
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	#pragma omp for
    	for (int j = 0; j < NX; j++) {
                for (int k = 0; k < NY; k++){ 
			float x = (j + 0.5) * DX;
			float y = (k + 0.5) * DY;
			int index2= j*NY+k;
        		T[index2] = 0.0;
        		if ((x > 0.1) && (x < 0.2) && (y < 0.2) && (y > 0.1)) {
            			T[index2] = 1.0;
        		}
		}
    	}

    /* Compute the analytical solution first */
	#pragma omp for
	for (int j = 0; j < NX; j++) {
                for (int k = 0; k < NY; k++){ 
                        float x = (j + 0.5) * DX;
                        float y = (k + 0.5) * DY;
                        int index2= j*NY+k;
                        A[index2] = 0.0;
                        if ((x>(0.1+U*T_FINAL)) && (x < 0.2+U*T_FINAL) && (y < 0.2+V*T_FINAL) && (y > 0.1+V*T_FINAL)) {
                                A[index2] = 1.0;
                        }
        	}
	}



    float time = 0.0;
    for (int timestep = 0; timestep < MAX_TIMESTEPS; timestep++) {
    
        // Compute fluxes
        Compute_Fluxes( T, F, W,D,time);

        // Update T using Fluxes F
//        printf("Computing state at time %g (using) timestep %g (after %d time steps)\n", time, DT, timestep);

        Update_State(F,T,W,D,Tnew);

        time = time + DT;
        if (time > T_FINAL) {
//            printf("Thread %d arrived at target time; stopping.\n", tid);
            break;
        }// else {
           // printf("Thread %d ran out of timesteps before reaching target time.\n", tid);
      //  }*/

    }

//    float Total_error = 0.0;
    #pragma omp single
    {
    for (int j = 0; j < NX; j++) {
                for (int k = 0; k < NY; k++){ 
			int index = j*NY+k;
			Total_error = Total_error + ((T[index] - A[index])*(T[index] - A[index]));
		} 
    }
    }

    // Find the average
    #pragma omp single
    {
    	Total_error = (double)(Total_error / N);
    	printf("Total error %g\n", Total_error);
	fprintf(pFile, "%d\t%g\n",N,Total_error); 
    }
    }//end of parallel
    fclose(pFile);
    if (DEBUG) printf("Saving T results\n");
    pFile = fopen("resultsT.txt", "w");
    for (int i = 0; i < NX; i++) {
    	for(int j = 0;j <NY;j++){
        	float X = (i+0.5)*DX;
                float Y = (j+0.5)*DY;
                int index = i*NY +j;
                fprintf(pFile, "%g\t%g\t%g\n",X, Y, T[index]);
        }
    }
    fclose(pFile);



    /* write fluxes to file: one line per interface (index, flux) */
    // Write U to file now

    /* cleanup */
    free(T);
    free(Tnew);
    free(F);
    free(A);
    free(W);
    free(D);
    return 0;
}
