#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 200          /* number of cells */
#define NIF (N+1)      /* number of interfaces (201) */

void Compute_Fluxes(int n, const double *u, double *F, double alpha)
{
    /* compute flux at all interfaces j = 0..n (n+1 interfaces)
       interface j is between left = (j-1) and right = j (mod n) */
    for (int j = 0; j <= n; ++j) {
        int left  = (j - 1 + n) % n;
        int right = j % n;
        double u_interface = (alpha > 0.0) ? u[left] : u[right];
        F[j] = alpha * u_interface;
    }
}

int main(void)
{
    double *u = malloc(sizeof(double) * N);
    double *F = malloc(sizeof(double) * NIF);
    if (!u || !F) {
        fprintf(stderr, "allocation failed\n");
        return 1;
    }

    /* initial condition: a smooth hump (Gaussian) in the domain [0,1) */
    for (int i = 0; i < N; ++i) {
        double x = (i + 0.5) / (double)N; /* cell center */
        u[i] = exp(-200.0 * (x - 0.5) * (x - 0.5));
    }

    double alpha = 0.1; /* given advection speed */
    Compute_Fluxes(N, u, F, alpha);

    /* write fluxes to file: one line per interface (index, flux) */
    FILE *fp = fopen("fluxes.dat", "w");
    if (!fp) {
        perror("fopen");
        free(u);
        free(F);
        return 1;
    }
    for (int j = 0; j < NIF; ++j) {
        fprintf(fp, "%d %.15e\n", j, F[j]);
    }
    fclose(fp);

    /* cleanup */
    free(u);
    free(F);

    return 0;
}