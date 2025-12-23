#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 200          /* number of cells */
#define NIF (N+1)      /* number of interfaces */
#define ALPHA 0.1    /* advection speed */
#define L 0.1          /* domain length */
#define DX (L / N)    /* cell size */
#define CFL 0.0000001
#define DT (CFL*DX / fabs(ALPHA)) /* time step size */   
#define T_FINAL 0.5     /* final time */
#define MAX_TIMESTEPS 500000


void Compute_Fluxes(int nif, const double *u, double *F, double alpha)
{
    /* compute flux at all interfaces j = 0..n (n+1 interfaces)
       interface j is between left = (j-1) and right = j (mod n) */

    // Only compute the inside interfaces; leave the outside 2 interfaces for later.   
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
    // Set dF/dx = 0 on left and right ends of our domain
    F[0] = F[1];
    F[nif-1] = F[nif-2];
}


void Update_State(int n, const double *F, double *u) {

    /*
    Solving
        du/dt + dF/dx = 0
    So
        U* = U - dt*dF/dx
    */
    for (int cell = 0; cell < n; cell++) {
        u[cell] = u[cell] - (DT/DX)*(F[cell+1] - F[cell]);
    }

}



int main(void)
{
    double *u = malloc(sizeof(double) * N);
    double *F = malloc(sizeof(double) * NIF);
    double *A = malloc(sizeof(double) * N);

    if (!u || !F) {
        fprintf(stderr, "allocation failed\n");
        return 1;
    }

    /* initial condition */
    for (int i = 0; i < N; ++i) {
        double x = (i + 0.5) * DX;
        u[i] = 0.0;
        if ((x > 0.02) && (x < 0.03)) {
            u[i] = 1.0;
        }
    }

    /* Compute the analytical solution first */

    for (int i = 0; i < N; ++i) {
        double x = (i + 0.5) * DX;
        A[i] = 0.0;
        if ((x > (0.02+ALPHA*T_FINAL)) && (x < 0.03+ALPHA*T_FINAL)) {
            A[i] = 1.0;
        }
    }



    double time = 0.0;
    for (int timestep = 0; timestep < MAX_TIMESTEPS; timestep++) {
    
        // Compute fluxes
        Compute_Fluxes(NIF, u, F, ALPHA);

        // Update U using Fluxes F
        printf("Computing state at time %g (using) timestep %g (after %d time steps)\n", time, DT, timestep);

        Update_State(N, F, u);

        time = time + DT;
        if (time > T_FINAL) {
            printf("Arrived at target time; stopping.\n");
            break;
        } else {
            printf("Ran out of timesteps before reaching target time.\n");
        }

    }

    double Total_error = 0.0;
    for (int i = 0; i < N; i++)
    {
        Total_error = Total_error + (u[i] - A[i])*(u[i] - A[i]);
    }
    // Find the average
    Total_error = (float)(Total_error / N);
    printf("Total error %g\n", Total_error); 

    





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