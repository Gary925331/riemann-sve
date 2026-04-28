#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NX 100          /* number of X cells */
#define NY 100
#define N (NX*NY)
#define NIF_X (NX+1)      /* number of interfaces */
#define NIF_Y (NY+1)
#define U 1    /* advection speed */
#define L 1          /* domain length */
#define H 1
#define DX L/NX    /* cell size */
#define DY H/NY
#define CFL DT*U/(DX*DX)
#define DT 0.0000001 /* time step size */   
#define T_FINAL 0.5     /* final time */
#define MAX_TIMESTEPS 5000000
#define NP 2
#define alpha 0.0005 //diffusion speed
#define DEBUG 1

void Compute_Fluxes(const float *T, float *F,float *W)
{
    /* compute flux at all interfaces j = 0..n (n+1 interfaces)
       interface j is between left = (j-1) and right = j (mod n) */

    // Only compute the inside interfaces; leave the outside 2 interfaces for later.
//	#pragma omp parallel
//	{
	#pragma omp for   
    	for (int j = 0; j < NIF_X; j++) {
		for (int k = 0; k < NIF_Y; k++){
			int index1 = j*NIF_Y+k;
			int index2 = k*NIF_X+j; 
			float Left_F = U*T[index2];
			float Right_F = U*T[index2+1];
			float Left,Right,Top,Bottom;
        		// Upwind scheme
        		if (U > 0) {
            			F[index2] = Left_F;
        		} else {
            			F[index2] = Right_F;
        		}
			if(j ==0){
                                Left = T[index2];
                        }else{
                                Left = T[index2];
                        }
			//Right
                        if (j == (NIF_X-1)) {
                                Right = T[index2];
                        }else{ 
                                Right = T[index2+1];
                        }
                                    
			//Bottom
                        if(k == 0){
                                Bottom = T[index1];
                        }else{
                                Bottom = T[index1];
                        }
                                   
			//Top
                        if(k == (NIF_Y-1)){
                                Top = T[index1];
			}else{
				Top = T[index1+1];
                        }

        		F[index2] += -alpha*((Right - Left) / DX);
			W[index1] = -alpha*((Top - Bottom) / DY);
        	}

       		//central flux
       		// F[j] = 0.5 * (left_F + right_F);


       		// Rusanov flux
       		// F[j] = 0.5 * (left_F + right_F) - 0.25*(right_u - left_u);


    	}
//	}//end of parallel
    	// Set dF/dx = 0 on left and right ends of our domain
    	
}


void Update_State(const float *F, float *T,float *W,float *Tnew) {

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
                        int index1 = j*NY+k;
                        int index2 = k*NX+j;
        		Tnew[index2] = T[index2] - (DT/DX)*(F[index2+1] - F[index2])-(DT/DY)*(W[index1+1]-W[index1]);
    		}
	}
	for (int j = 0; j < NX; j++) {
                for (int k = 0; k < NY; k++){
                        int index1 = j*NY+k;
                        int index2 = k*NX+j; 
                        T[index2] = Tnew[index2];
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
	T = (float*)malloc(NIF_X*NIF_Y*sizeof(float));
	Tnew = (float*)malloc(NIF_X*NIF_Y*sizeof(float));
	F = (float*)malloc(NIF_X*NIF_Y*sizeof(float));
	A = (float*)malloc(NIF_X*NIF_Y*sizeof(float));
	W = (float*)malloc(NIF_X*NIF_Y*sizeof(float));
    	if (!T || !F) {
        	fprintf(stderr, "allocation failed\n");
        	return 1;
    	}

    /* initial condition */
	omp_set_num_threads(NP);
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
        		if ((x > 0.2) && (x < 0.4) && (y < 0.6) && (y > 0.4)) {
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
                        if ((0.2+U*T_FINAL) && (x < 0.4+U*T_FINAL) && (y < 0.6) && (y > 0.4)) {
                                A[index2] = 1.0;
                        }
        	}
	}



    float time = 0.0;
    for (int timestep = 0; timestep < MAX_TIMESTEPS; timestep++) {
    
        // Compute fluxes
        Compute_Fluxes( T, F, W);

        // Update T using Fluxes F
//        printf("Computing state at time %g (using) timestep %g (after %d time steps)\n", time, DT, timestep);

        Update_State(F,T,W,Tnew);

        time = time + DT;
/*        if (time > T_FINAL) {
            printf("Thread %d arrived at target time; stopping.\n", tid);
            break;
        } else {
            printf("Thread %d ran out of timesteps before reaching target time.\n", tid);
        }*/

    }

    float Total_error = 0.0;
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
    	Total_error = (float)(Total_error / N);
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
                fprintf(pFile, "%g\t%g\t%g\n", X, Y, T[index]);
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
    return 0;
}
