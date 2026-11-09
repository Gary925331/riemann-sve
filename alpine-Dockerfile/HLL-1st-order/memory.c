#include <stdlib.h>
#include <stdio.h>

void Allocate_memory(float **u,float **v,float **mass_F,float **momentum_F_X,float **momentum_F_Y,float **mass_G,float **momentum_G_X,float **momentum_G_Y,float **mass,
float **momentum_X,float **momentum_Y,float **h,int N,int NIF){
	size_t alignment = 32;
	int err = 0;
	err |= posix_memalign((void**)u,            alignment, N*sizeof(float));
	err |= posix_memalign((void**)v,            alignment, N*sizeof(float));
	err |= posix_memalign((void**)mass_F,       alignment, NIF*sizeof(float));
	err |= posix_memalign((void**)momentum_F_X, alignment, NIF*sizeof(float));
	err |= posix_memalign((void**)momentum_F_Y, alignment, NIF*sizeof(float));
	err |= posix_memalign((void**)mass_G,       alignment, NIF*sizeof(float));
        err |= posix_memalign((void**)momentum_G_X, alignment, NIF*sizeof(float));
	err |= posix_memalign((void**)momentum_G_Y, alignment, NIF*sizeof(float));
	err |= posix_memalign((void**)mass,         alignment, NIF*sizeof(float));
	err |= posix_memalign((void**)momentum_X,   alignment, NIF*sizeof(float));
	err |= posix_memalign((void**)momentum_Y,   alignment, NIF*sizeof(float));
	err |= posix_memalign((void**)h,            alignment, N*sizeof(float));
	if (err != 0) {
		fprintf(stderr, "posix_memalign 配置記憶體失敗\n");
		exit(EXIT_FAILURE);
	}
	/*h_slope_X = (float*)malloc(N*sizeof(float));
	*u_slope_X_X = (float*)malloc(N*sizeof(float));
	*v_slope_X_Y = (float*)malloc(N*sizeof(float));
	*h_slope_Y = (float*)malloc(N*sizeof(float));
        *u_slope_Y_X = (float*)malloc(N*sizeof(float));
	*v_slope_Y_Y = (float*)malloc(N*sizeof(float));*/
}
void Free_memory(float **u,float **v,float **mass_F,float **momentum_F_X,float **momentum_F_Y,float **mass_G,float **momentum_G_X,float **momentum_G_Y,float **mass,float **momentum_X,float **momentum_Y,float **h){
	free(*u);
	free(*v);
	free(*mass_F);
	free(*momentum_F_X);
	free(*momentum_F_Y);
	free(*mass_G);
        free(*momentum_G_X);
	free(*momentum_G_Y);
	free(*mass);
	free(*momentum_X);
	free(*momentum_Y);
	free(*h);
	/*free(*h_slope_X);
	free(*u_slope_X_X);
	free(*v_slope_X_Y);
	free(*h_slope_Y);
        free(*u_slope_Y_X);
	free(*v_slope_Y_Y);*/
}
