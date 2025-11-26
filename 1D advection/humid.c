#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// --- Simulation Parameters for Humid Air in Pipe ---
#define N_CELLS 200          // Number of cells (N)
#define N_FLUXES (N_CELLS + 1) // Number of interfaces (N+1)
#define ALPHA 0.1f           // Advection velocity (Flow speed: 0.1 m/s)
#define DOMAIN_LENGTH 0.1f   // Length of the pipe (0.1 m)
#define FINAL_TIME 0.5f      // Target simulation time (0.5 s)

// Target times for output
#define T_HALF 0.25f
#define T_FINAL FINAL_TIME

// --- CFL Condition Parameter ---
#define CFL_NUM 0.5f         // Chosen CFL number (must be <= 1.0 for Upwind scheme)

// Calculated parameters
#define DX (DOMAIN_LENGTH / N_CELLS) // Cell size (Delta x = 0.0005 m)
// DT calculated based on the explicitly defined CFL_NUM: DT = CFL_NUM * DX / |ALPHA|
#define DT (CFL_NUM * DX / fabsf(ALPHA)) 
// Total number of steps needed to reach FINAL_TIME
#define N_STEPS ((int)(FINAL_TIME / DT)) 

/**
 * @brief Saves the array data to a CSV file.
 * @param filename The name of the file to write.
 * @param data The array of float values (humidity).
 * @param size The number of elements in the array.
 * @param title A header string for the file (e.g., "Humidity_Values").
 */
void Save_To_File(const char *filename, const float data[], int size, const char *title) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    fprintf(file, "Position_x,Cell_Index,%s\n", title);
    for (int i = 0; i < size; i++) {
        // Calculate the cell center position for plotting
        float x_pos = (i + 0.5f) * DX; 
        fprintf(file, "%.6f,%d,%.10f\n", x_pos, i, data[i]);
    }

    fclose(file);
    printf("Successfully saved data to %s\n", filename);
}

/**
 * @brief Initializes the cell center array U (humidity) with a pulse between 0.02m and 0.03m.
 * @param U The array of cell values (humidity).
 * @param N The number of cells.
 */
void Initialize_U(float U[], int N) {
    // Determine the cell indices corresponding to the physical boundaries
    // Start index for 0.02m: I_start = floor(0.02 / DX)
    // End index for 0.03m: I_end = floor(0.03 / DX)
    const int I_START = (int)(0.02f / DX);
    const int I_END = (int)(0.03f / DX);

    // Initialize all cells to 0.0 (0% humidity)
    for (int i = 0; i < N; i++) {
        U[i] = 0.0f;
    }

    // Set the humidity pulse to 1.0 (100% humidity)
    // For DX = 0.0005: I_START = 40, I_END = 60. Pulse is in cells 40 through 59.
    for (int i = I_START; i < I_END; i++) {
        if (i < N) {
            U[i] = 1.0f;
        }
    }
    printf("Initial humidity pulse set from x=%.4f to x=%.4f (Cells %d to %d).\n", 
           I_START * DX, I_END * DX, I_START, I_END - 1);
}

/**
 * @brief Computes the numerical flux F at all N+1 cell interfaces using the Upwind scheme.
 * F[i] is the flux across the interface between cell i-1 and cell i.
 * Since alpha = 0.1 > 0, we use the value from the left (upstream) cell.
 */
void Compute_Fluxes(int N, const float U[], float F[], float alpha) {
    // 1. Calculate the interior and right boundary fluxes (F[1] to F[N])
    // Loop runs from i=1 to N. F[i] is between U[i-1] and U[i].
    for (int i = 1; i <= N; i++) {
        if (alpha > 0.0f) {
            // Flow to the right (alpha > 0): Upwind source is U[i-1] (to the left)
            F[i] = alpha * U[i - 1];
        } else {
            // Flow to the left (alpha < 0): Outflow boundary used here for completeness.
            float u_right = (i == N) ? U[N - 1] : U[i];
            F[i] = alpha * u_right;
        }
    }

    // 2. Calculate the left boundary flux (F[0])
    // F[0] is the flux entering cell 0 from the left.
    if (alpha > 0.0f) {
        // Inflow boundary: Assume zero inflow (U_BC = 0.0, 0% humidity entering)
        F[0] = 0.0f;
    } else {
        // Outflow boundary: Flux is determined by the upstream cell (U[0])
        F[0] = alpha * U[0];
    }
}

/**
 * @brief Updates the cell values U using the computed fluxes F (Forward Euler in time).
 * Governing equation: U_i^{n+1} = U_i^n - (DT/DX) * (F_{i+1/2} - F_{i-1/2})
 */
void Update_U(int N, float U[], const float F[], float dt, float dx) {
    float dt_over_dx = dt / dx;

    // Loop over all cells (i = 0 to N-1)
    for (int i = 0; i < N; i++) {
        // F[i] is the flux *entering* cell i (F_{i-1/2})
        // F[i+1] is the flux *leaving* cell i (F_{i+1/2})
        float flux_difference = F[i + 1] - F[i];
        
        // Update U_i
        U[i] = U[i] - dt_over_dx * flux_difference;
    }
}

int main() {
    printf("Starting 1D Humid Air Advection Solver (Upwind Scheme)\n");
    printf("Pipe Length (L): %.2f m\n", DOMAIN_LENGTH);
    printf("Flow Speed (ALPHA): %.2f m/s\n", ALPHA);
    printf("CFL Number: %.2f\n", CFL_NUM);
    printf("Total Simulation Time: %.2f s\n", FINAL_TIME);
    printf("DX: %.6f m, DT: %.6f s, Total Steps: %d\n\n", DX, DT, N_STEPS);

    // 1. Allocate memory for cell values (U) and fluxes (F)
    float *U = (float *)malloc(N_CELLS * sizeof(float));
    float *F = (float *)malloc(N_FLUXES * sizeof(float));

    if (U == NULL || F == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        if (U) free(U);
        if (F) free(F);
        return 1;
    }

    // 2. Initialize the humidity profile
    Initialize_U(U, N_CELLS);
    float current_time = 0.0f;
    
    // Save initial state for comparison (t=0s)
    Save_To_File("humidity_t_0_00s.csv", U, N_CELLS, "Humidity_t_0.00s");
    printf("Profile saved at t=0.00s.\n");

    // 3. Main Time Stepping Loop
    for (int step = 0; step < N_STEPS; step++) {
        // a. Compute the fluxes F based on the current U
        Compute_Fluxes(N_CELLS, U, F, ALPHA);

        // b. Update the cell values U based on the fluxes
        Update_U(N_CELLS, U, F, DT, DX);

        current_time += DT;

        // Save profile at t=0.25s
        if (current_time >= T_HALF - DT/2.0f && current_time < T_HALF + DT/2.0f) {
             // Create a copy to save the state precisely at T_HALF
            float *U_half = (float *)malloc(N_CELLS * sizeof(float));
            if (U_half) {
                for (int i = 0; i < N_CELLS; i++) U_half[i] = U[i];
                Save_To_File("humidity_t_0_25s.csv", U_half, N_CELLS, "Humidity_t_0.25s");
                free(U_half);
                printf("Profile saved at t=0.25s (Step %d).\n", step + 1);
            }
        }
    }

    // 4. Save the final computed U array to a file (t=0.5s)
    // Note: Due to floating point math, current_time might not be exactly 0.5, but close.
    Save_To_File("humidity_t_0_50s.csv", U, N_CELLS, "Humidity_t_0.50s");
    printf("Profile saved at t=0.50s.\n");

    // 5. Clean up
    free(U);
    free(F);
    
    printf("\nSimulation finished. Three time profiles saved.\n");

    return 0;
}