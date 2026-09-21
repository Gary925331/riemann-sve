#include <stdlib.h>
#include <stdio.h>

void Allocate_memory(float **u,float **v,float **s,float **mass_F,float **momentum_F_X,float **momentum_F_Y,float **mass_G,float **momentum_G_X,float **momentum_G_Y,float **mass,float **momentum_X,float **momentum_Y,float **h,float **h_slope_X,float **u_slope_X_X,float **v_slope_X_Y,float **h_slope_Y,float **u_slope_Y_X,float **v_slope_Y_Y,int N,int NIF){
	*u = (float*)malloc(N*sizeof(float));
	*v = (float*)malloc(N*sizeof(float));
	*s = (float*)malloc(N*sizeof(float));
	*mass_F = (float*)malloc(NIF*sizeof(float));
	*momentum_F_X = (float*)malloc(NIF*sizeof(float));
	*momentum_F_Y = (float*)malloc(NIF*sizeof(float));
	*mass_G = (float*)malloc(NIF*sizeof(float));
        *momentum_G_X = (float*)malloc(NIF*sizeof(float));
	*momentum_G_Y = (float*)malloc(NIF*sizeof(float));
	*mass = (float*)malloc(NIF*sizeof(float));
	*momentum_X = (float*)malloc(NIF*sizeof(float));
	*momentum_Y = (float*)malloc(NIF*sizeof(float));
	*h = (float*)malloc(N*sizeof(float));
	*h_slope_X = (float*)malloc(N*sizeof(float));
	*u_slope_X_X = (float*)malloc(N*sizeof(float));
	*v_slope_X_Y = (float*)malloc(N*sizeof(float));
	*h_slope_Y = (float*)malloc(N*sizeof(float));
        *u_slope_Y_X = (float*)malloc(N*sizeof(float));
	*v_slope_Y_Y = (float*)malloc(N*sizeof(float));
}
void Free_memory(float **u,float **v,float **s,float **mass_F,float **momentum_F_X,float **momentum_F_Y,float **mass_G,float **momentum_G_X,float **momentum_G_Y,float **mass,float **momentum_X,float **momentum_Y,float **h,float **h_slope_X,float **u_slope_X_X,float **v_slope_X_Y,float **h_slope_Y,float **u_slope_Y_X,float **v_slope_Y_Y){
	free(*u);
	free(*v);
	free(*s);
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
	free(*h_slope_X);
	free(*u_slope_X_X);
	free(*v_slope_X_Y);
	free(*h_slope_Y);
        free(*u_slope_Y_X);
	free(*v_slope_Y_Y);
}
