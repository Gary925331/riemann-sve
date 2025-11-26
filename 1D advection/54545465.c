#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// --- Simulation Parameters ---
#define N_CELLS 200          // Number of cells (N)
#define N_FLUXES (N_CELLS + 1) // Number of interfaces (N+1)
#define ALPHA 0.1f           // Advection velocity (alpha)
#define DOMAIN_LENGTH 1.0f   // Length of the domain
#define N_STEPS 100          // Number of time steps to run

// Calculated parameters
#define DX (DOMAIN_LENGTH / N_CELLS) // Cell size (Delta x)
// DT calculated for CFL = 0.5: DT = CFL * DX / |ALPHA|
#define DT (0.5f * DX / fabsf(ALPHA)) 

/**
 * @brief Saves the array data to a CSV file.
 * @param filename The name of the file to write.
 * @param data The array of float values.
 * @param size The number of elements in the array.
 * @param title A header string for the file (e.g., "U_Values").
 */
void Save_To_File(const char *filename, const float data[], int size, const char *title) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    fprintf(file, "Cell_Index,%s\n", title);
    for (int i = 0; i < size; i++) {
        fprintf(file, "%d,%.10f\n", i, data[i]);
    }

    fclose(file);
    printf("Successfully saved data to %s\n", filename);
}

/**
 * @brief Initializes the cell center array U with a simple pulse.
 * @param U The array of cell values.
 * @param N The number of cells.
 */
void Initialize_U(float U[], int N) {
    // Initialize all cells to 0.0
    for (int i = 0; i < N; i++) {
        U[i] = 0.0f;
    }

    // Set a simple square pulse in the center of the domain (cells 50 to 100)
    for (int i = 50; i < 100; i++) {
        U[i] = 1.0f;
    }
}

/**
 * @brief Computes the numerical flux F at all N+1 cell interfaces using the Upwind scheme.
 * F[i] is the flux across the interface between cell i-1 and cell i.
 * * Upwind Rule:
 * If alpha > 0 (flow to the right), F[i] = alpha * U[i-1] (flux from the left).
 * If alpha < 0 (flow to the left), F[i] = alpha * U[i] (flux from the right).
 *
 * @param N The number of cells (size of U).
 * @param U The cell-centered value array (size N).
 * @param F The interface flux array (size N+1).
 * @param alpha The constant advection velocity.
 */
void Compute_Fluxes(int N, const float U[], float F[], float alpha) {
    // 1. Calculate the interior and right boundary fluxes (F[1] to F[N])
    // Loop runs from i=1 to N. F[i] is between U[i-1] and U[i].
    for (int i = 1; i <= N; i++) {
        if (alpha > 0.0f) {
            // Flow to the right (alpha > 0): Upwind source is U[i-1] (to the left)
            F[i] = alpha * U[i - 1];
        } else {
            // Flow to the left (alpha < 0): Upwind source is U[i] (to the right)
            // Note: For i=N, U[N] would be required (ghost cell).
            // We implement an outflow boundary condition here: U[N] = U[N-1]
            float u_right = (i == N) ? U[N - 1] : U[i];
            F[i] = alpha * u_right;
        }
    }

    // 2. Calculate the left boundary flux (F[0])
    // F[0] is the flux entering cell 0 from the left.
    if (alpha > 0.0f) {
        // Inflow boundary: Assume zero inflow (U_BC = 0.0)
        F[0] = 0.0f;
    } else {
        // Outflow boundary: Flux is determined by the upstream cell (U[0])
        F[0] = alpha * U[0];
    }
}

/**
 * @brief Updates the cell values U using the computed fluxes F (Forward Euler in time).
 * Governing equation: U_i^{n+1} = U_i^n - (DT/DX) * (F_{i+1/2} - F_{i-1/2})
 *
 * @param N The number of cells.
 * @param U The cell-centered value array (will be modified).
 * @param F The interface flux array.
 * @param dt Time step size.
 * @param dx Cell size.
 */
void Update_U(int N, float U[], const float F[], float dt, float dx) {
    float dt_over_dx = dt / dx;

    // Loop over all internal cells (i = 0 to N-1)
    for (int i = 0; i < N; i++) {
        // The flux *entering* cell i is F[i] (F_{i-1/2})
        // The flux *leaving* cell i is F[i+1] (F_{i+1/2})
        float flux_difference = F[i + 1] - F[i];
        
        // Update U_i
        U[i] = U[i] - dt_over_dx * flux_difference;
    }
}

int main() {
    printf("Starting 1D Advection Solver (Upwind Scheme)\n");
    printf("N=%d, DX=%.4f, DT=%.4f, ALPHA=%.1f\n", N_CELLS, DX, DT, ALPHA);

    // 1. Allocate memory for cell values and fluxes (using float)
    float *U = (float *)malloc(N_CELLS * sizeof(float));
    float *F = (float *)malloc(N_FLUXES * sizeof(float));

    if (U == NULL || F == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        // Ensure allocated memory is freed before exiting
        if (U) free(U);
        if (F) free(F);
        return 1;
    }

    // 2. Initialize the cell values
    Initialize_U(U, N_CELLS);
    
    // Save initial state for comparison
    Save_To_File("initial_u_profile.csv", U, N_CELLS, "Initial_U");

    // 3. Main Time Stepping Loop
    for (int step = 0; step < N_STEPS; step++) {
        // a. Compute the fluxes F based on the current U
        Compute_Fluxes(N_CELLS, U, F, ALPHA);

        // b. Update the cell values U based on the fluxes
        Update_U(N_CELLS, U, F, DT, DX);

        if ((step + 1) % (N_STEPS / 10) == 0) {
            printf("Step %d completed.\n", step + 1);
        }
    }

    // 4. Save the final computed U array to a file
    Save_To_File("final_u_profile.csv", U, N_CELLS, "Final_U");

    // 5. Clean up
    free(U);
    free(F);
    
    printf("Simulation finished and final profile saved.\n");

    return 0;
}