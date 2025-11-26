#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// --- Constants and Definitions ---
#define N_CELLS 200      // Number of cells (N)
#define N_INTERFACES 201 // Number of interfaces/fluxes (N+1)
#define ALPHA 0.1f       // Advection velocity (alpha), using 'f' suffix for float

/**
 * @brief Computes the numerical flux F at all cell interfaces using the Upwind Scheme.
 * * This version uses float for all calculations.
 *
 * @param U Array of cell-centered values (size N_CELLS).
 * @param F Array to store the interface fluxes (size N_INTERFACES).
 * @param alpha The advection velocity (float).
 */
void Compute_Fluxes(float *U, float *F, float alpha) {
    int i_left, i_right;
    
    // Calculate fluxes for all N+1 interfaces (j = 0 to N_CELLS)
    for (int j = 0; j <= N_CELLS; j++) {
        
        // The interface j sits between cell i_left and cell i_right.
        
        // Periodic Boundary Conditions:
        if (j == 0) {
            // Left boundary flux (F[0]): Upstream is the last cell U[N_CELLS-1]
            i_left = N_CELLS - 1;
            i_right = 0;
        } else if (j == N_CELLS) {
            // Right boundary flux (F[N]): Upstream of F[N] is U[N_CELLS-1].
            i_left = N_CELLS - 1;
            i_right = 0; // The "right" cell wraps back to cell 0
        } else {
            // Interior interfaces:
            i_left = j - 1;
            i_right = j;
        }
        
        // --- Upwind Scheme Implementation ---
        if (alpha > 0.0f) { // Use 'f' suffix for float comparison
            // Flow is from left to right. Flux uses the value from the cell on the left.
            F[j] = alpha * U[i_left];
        } else if (alpha < 0.0f) { // Use 'f' suffix for float comparison
            // Flow is from right to left. Flux uses the value from the cell on the right.
            F[j] = alpha * U[i_right];
        } else {
            // alpha == 0, no advection.
            F[j] = 0.0f;
        }
    }
}

/**
 * @brief Initializes the cell values U with a simple hat function (float).
 * @param U Array of cell-centered values.
 */
void Initialize_U(float *U) {
    int start = N_CELLS / 4;
    int end = N_CELLS / 2;
    
    for (int i = 0; i < N_CELLS; i++) {
        if (i >= start && i < end) {
            // Simple hat/block function for visualization
            U[i] = 1.0f; // Use 'f' suffix
        } else {
            U[i] = 0.0f; // Use 'f' suffix
        }
    }
}

/**
 * @brief Saves the flux array F to a text file (float).
 * @param F Array of interface fluxes.
 * @param filename Name of the file to save.
 */
void Save_Fluxes(float *F, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file for writing");
        return;
    }
    
    // Save N_INTERFACES (201) flux values
    for (int j = 0; j < N_INTERFACES; j++) {
        // Using %f for float output
        fprintf(file, "%d, %.8f\n", j, F[j]);
    }
    
    fclose(file);
    printf("Successfully saved %d fluxes to %s\n", N_INTERFACES, filename);
}

// --- Main Function ---
int main() {
    // 1. Array Declarations and Dynamic Allocation
    // U: N_CELLS (200) cell-centered values
    float *U = (float *)malloc(N_CELLS * sizeof(float));
    // F: N_INTERFACES (201) interface flux values
    float *F = (float *)malloc(N_INTERFACES * sizeof(float));

    // Check for memory allocation failure
    if (U == NULL || F == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        // Ensure allocated memory is freed before exiting if one allocation succeeded
        if (U) free(U);
        if (F) free(F);
        return 1; // Return non-zero to indicate error
    }
    
    // 2. Initialize cell values (U)
    Initialize_U(U);
    printf("Initialized %d cell values (U).\n", N_CELLS);

    // 3. Compute Fluxes 
    printf("Computing fluxes F with alpha = %.2f...\n", ALPHA);
    Compute_Fluxes(U, F, ALPHA);

    // 4. Save results
    Save_Fluxes(F, "flux_output.csv");
    
    // Optional: Print a few values for verification
    printf("\nFirst few fluxes (F):\n");
    for (int j = 0; j < 5; j++) {
        printf("F[%d] = %.8f\n", j, F[j]);
    }
    printf("...\n");
    printf("Last few fluxes (F):\n");
    for (int j = N_INTERFACES - 5; j < N_INTERFACES; j++) {
        printf("F[%d] = %.8f\n", j, F[j]);
    }

    // 5. Free allocated memory
    free(U);
    free(F);

    return 0;
}