#include <stdio.h>
#include <stdlib.h> // Required for malloc, free, exit

#define NO_CELLS 200
#define NO_INTERFACES (NO_CELLS + 1)
#define NO_TIME_STEPS 100000

// --- Helper Function ---
// Returns a pointer to the new memory block
float* allocate_memory(int number_of_elements) {
    float* ptr = (float*)malloc(number_of_elements * sizeof(float));
    
    if (ptr == NULL) {
        printf("Error: Memory allocation failed.\n");
        exit(1);
    }
    
    return ptr;
}

int main() {
    // 1. Parameter Definitions
    float L = 1.0f;
    float U = 0.5f;
    float ALPHA = 0.1f;
    float DX = L / NO_CELLS;
    float CFL = 0.0001f;
    float DT = CFL * DX / ALPHA;

    // 2. Allocate Memory (Simple assignment)
    float *H = allocate_memory(NO_CELLS);
    float *H_INIT = allocate_memory(NO_CELLS);
    float *H_new = allocate_memory(NO_CELLS);
    float *F = allocate_memory(NO_INTERFACES);

    // 3. Initialize Arrays to 0.0 (Using loops, not memset)
    for (int i = 0; i < NO_CELLS; i++) {
        H[i] = 0.0f;
        H_INIT[i] = 0.0f;
        H_new[i] = 0.0f;
    }
    for (int i = 0; i < NO_INTERFACES; i++) {
        F[i] = 0.0f;
    }

    // 4. Setup Initial Pulse
    int start = (int)(0.2f / DX);
    int end = (int)(0.4f / DX);

    for (int i = start; i < end; i++) {
        H[i] = 1.0f;
    }

    printf("Initial pulse in cells %d to %d\n", start, end - 1);

    // Copy initial state
    for (int i = 0; i < NO_CELLS; i++) {
        H_INIT[i] = H[i];
    }

    // 5. Time Stepping Loop
    for (int step = 0; step < NO_TIME_STEPS; step++) {
        
        if (step % 10000 == 0) {
            printf("Computing step %d\n", step);
        }

        // Flux Calculation
        for (int i = 1; i < NO_INTERFACES - 1; i++) {
            F[i] = 0.0f;

            if (U > 0) {
                F[i] += U * H[i - 1];
            } else {
                F[i] += U * H[i];
            }
            
            F[i] += -ALPHA * (H[i] - H[i - 1]) / DX;
        }

        // Update Cell Averages
        for (int i = 0; i < NO_CELLS; i++) {
            H_new[i] = H[i] - (DT / DX) * (F[i + 1] - F[i]);
        }

        // Update H for next step
        for (int i = 0; i < NO_CELLS; i++) {
            H[i] = H_new[i];
        }
    }

    // 6. Output Result
    FILE *fp = fopen("results.csv", "w");
    if (fp != NULL) {
        fprintf(fp, "x,H\n");
        for (int i = 0; i < NO_CELLS; i++) {
            fprintf(fp, "%f,%f\n", i * DX, H[i]);
        }
        fclose(fp);
        printf("Data saved to results.csv\n");
    }

    // 7. Clean up
    free(H);
    free(H_INIT);
    free(H_new);
    free(F);

    return 0;
}