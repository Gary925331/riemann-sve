#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define N 200          /* number of cells */
#define NIF (N+1)      /* number of interfaces */
#define ALPHA 0.1    /* advection speed */
#define L 0.1          /* domain length */
#define DX (L / N)    /* cell size */
#define CFL 0.0000001
#define DT (CFL*DX / fabs(ALPHA)) /* time step size */   
#define T_FINAL 0.0005     /* final time */
#define MAX_TIMESTEPS 5000000
#define NP 2


void Compute_Fluxes(int nif, const float *u, float *F,float alpha)
{
    /* compute flux at all interfaces j = 0..n (n+1 interfaces)
       interface j is between left = (j-1) and right = j (mod n) */

    // Only compute the inside interfaces; leave the outside 2 interfaces for later.
//	#pragma omp parallel
//	{
	#pragma omp for   
    	for (int j = 1; j < (nif-1); j++) {
        	int left  = j - 1;
        	int right = j;
        
        	double left_F = ALPHA*u[left];
        	double right_F = ALPHA*u[right];
        	double left_u = u[left];
        	double right_u = u[right];
        
        	// Upwind scheme
        
        	if (ALPHA > 0) {
            		F[j] = left_F;
        	} else {
            		F[j] = right_F;
        	}

        	F[j] += -ALPHA*(u[j] - u[j-1]) / DX;
        

       		//central flux
       		// F[j] = 0.5 * (left_F + right_F);


       		// Rusanov flux
       		// F[j] = 0.5 * (left_F + right_F) - 0.25*(right_u - left_u);


    	}
//	}//end of parallel
    	// Set dF/dx = 0 on left and right ends of our domain
    	F[0] = F[1];
    	F[nif-1] = F[nif-2];
}


void Update_State(int n, const float *F, float *u) {

    /*
    Solving
        du/dt + dF/dx = 0
    So
        U* = U - dt*dF/dx
    */
//	#pragma omp parallel
//	{
	#pragma omp for
    	for (int cell = 0; cell < n; cell++) {
        	u[cell] = u[cell] - (DT/DX)*(F[cell+1] - F[cell]);
    	}
//	}

}



int main(void)
{
	float *u;
	float *F;
	float *A;
	u = (float*)malloc(N*sizeof(float));
	F = (float*)malloc(NIF*sizeof(float));
	A = (float*)malloc(N*sizeof(float));
    	if (!u || !F) {
        	fprintf(stderr, "allocation failed\n");
        	return 1;
    	}

    /* initial condition */
	omp_set_num_threads(NP);
	#pragma omp parallel
	{
	int tid = omp_get_thread_num();
	#pragma omp for
    	for (int i = 0; i < N; ++i) {
        	double x = (i + 0.5) * DX;
        	u[i] = 0.0;
        	if ((x > 0.02) && (x < 0.03)) {
            	u[i] = 1.0;
        	}
    	}

    /* Compute the analytical solution first */
	#pragma omp for
    	for (int i = 0; i < N; ++i) {
        	double x = (i + 0.5) * DX;
        	A[i] = 0.0;
        	if ((x > (0.02+ALPHA*T_FINAL)) && (x < 0.03+ALPHA*T_FINAL)) {
            		A[i] = 1.0;
        	}
    	}
//	}//end of parallel



    float time = 0.0;
    for (int timestep = 0; timestep < MAX_TIMESTEPS; timestep++) {
    
        // Compute fluxes
        Compute_Fluxes(NIF, u, F, ALPHA);

        // Update U using Fluxes F
        printf("Computing state at time %g (using) timestep %g (after %d time steps)\n", time, DT, timestep);

        Update_State(N, F, u);

        time = time + DT;
        if (time > T_FINAL) {
            printf("Thread %d arrived at target time; stopping.\n", tid);
            break;
        } else {
            printf("Thread %d ran out of timesteps before reaching target time.\n", tid);
        }

    }

    float Total_error = 0.0;
    #pragma omp single
    for (int i = 0; i < N; i++)
    {
        Total_error = Total_error + (u[i] - A[i])*(u[i] - A[i]);
    }
    // Find the average
    #pragma omp single
    {
    	Total_error = (float)(Total_error / N);
    	printf("Total error %g\n", Total_error); 
    }
    }//end of parallel





    /* write fluxes to file: one line per interface (index, flux) */
    FILE *fp = fopen("fluxes.dat", "w");
    for (int j = 0; j < NIF; ++j) {
        fprintf(fp, "%d %.15e\n", j, F[j]);
    }
    fclose(fp);

    // Write U to file now
    fp = fopen("results.dat", "w");
    for (int cell = 0; cell < N; cell++) {
        double x = (cell+0.5)*DX;
        fprintf(fp, "%g\t%g\n", x, u[cell]);
    }
    fclose(fp);

    /* cleanup */
    free(u);
    free(F);
    free(A);

    return 0;
}
