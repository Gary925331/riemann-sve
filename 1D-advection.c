#include <stdio.h>
#include <stdlib.h>

int main() {
    const int N = 200;// 陣列大小
    float *U, *x; // 宣告兩個浮點陣列指標
    // 配置記憶體
    U = (float*)malloc(N * sizeof(float));
    x = (float*)malloc(N * sizeof(float));
    // 確認是否配置成功
    if (U == NULL || x == NULL) {
        printf("Memory allocation failed. Program aborted.\n");
        if (U) free(U);
        if (x) free(x);
        return 1;
    }
    // 初始化測試（非必要，可省略
    int i;   // <-- C89 必須先宣告變數
    for (i = 0; i < N; i++) {
        U[i] = 0.0f;
        x[i] = (float)i;
    }

    printf("Successfully allocated and initialized U and x arrays with %d elements.\n", N);
    // 釋放記憶體
    free(U);
    free(x);

    printf("Memory successfully freed.\n");
    return 0;
}
