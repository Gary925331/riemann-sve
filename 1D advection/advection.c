#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

    const int N = 200;
    float *U, *X;
    int i;

    // 配置記憶體
    U = (float*)malloc(N * sizeof(float));
    X = (float*)malloc(N * sizeof(float));

    if ((U == NULL) || (X == NULL)) {
        printf("Unable to allocate memory for U or X. Aborting\n");
        if (U) free(U);
        if (X) free(X);
        return 1;
    }

    // 初始化資料
    for (i = 0; i < N; i++) {
        X[i] = (float)i;
        U[i] = sqrtf(X[i] + 2 * X[i]);  // 等同於 sqrt(3*i)
        printf("Sqrt of (%g+%g) = %e\n", X[i], 2 * X[i], U[i]);
    }

    // 釋放記憶體
    free(U);
    free(X);

    return 0;
}
