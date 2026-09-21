#include <stdlib.h>
#include <math.h>
#include <omp.h>

void ghost(float *u,float *v,float *s,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,
float *momentum_Y,float *h,float *h_slope_X,float *u_slope_X_X,float *v_slope_X_Y,float *h_slope_Y,float *u_slope_Y_X,float *v_slope_Y_Y){	

		#pragma omp for
                for(int j = 0; j < NY+2; j++){
                        mass[0*(NY+2)+j] = mass[1*(NY+2)+j];
                        momentum_X[0*(NY+2)+j] = momentum_X[1*(NY+2)+j]; 
                        momentum_Y[0*(NY+2)+j] = momentum_Y[1*(NY+2)+j];

                        mass[(NX+1)*(NY+2)+j] = mass[NX*(NY+2)+j];
                        momentum_X[(NX+1)*(NY+2)+j] = momentum_X[NX*(NY+2)+j]; 
                        momentum_Y[(NX+1)*(NY+2)+j] = momentum_Y[NX*(NY+2)+j];

                }
                #pragma omp for
                for(int i = 0; i < NX+2; i++){
                        mass[i*(NY+2)+0] = mass[i*(NY+2)+1];
                        momentum_X[i*(NY+2)+0] = momentum_X[i*(NY+2)+1];
                        momentum_Y[i*(NY+2)+0] = -momentum_Y[i*(NY+2)+1]; // 下牆反彈

                        mass[i*(NY+2)+NY+1] = mass[i*(NY+2)+NY];
                        momentum_X[i*(NY+2)+NY+1] = momentum_X[i*(NY+2)+NY];
                        momentum_Y[i*(NY+2)+NY+1] = -momentum_Y[i*(NY+2)+NY]; // 上牆反彈
                }
                #pragma omp for
                for(int i = 0; i < NX+2; i++){
                        for (int j = 0; j < NY+2; j++){
                                int index = i*(NY+2)+j;
                                h[index] = mass[index];
                                u[index] = momentum_X[index]/mass[index];
                                v[index] = momentum_Y[index]/mass[index];
                        }
                }
}
