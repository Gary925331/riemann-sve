#include <stdlib.h>
#include <math.h>
#include <omp.h>
void state(float *u,float *v,float *s,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,
float *momentum_Y,float *h,float *h_slope_X,float *u_slope_X_X,float *v_slope_X_Y,float *h_slope_Y,float *u_slope_Y_X,float *v_slope_Y_Y,int NX,int NY,float DX,float DY,float DT){
		#pragma omp for
                for(int i = 1;i < NX+1;i++){
                        for (int j = 1; j < NY+1; j++){
                                int index = i*(NY+2)+j;
                                //int index1 = i*(NIF_Y+2)+j;
                                if (i == NX/2 +1 && (j <= NY*96/200 || j >= NY*171/200)) {
                                        continue; 
                                }
                                mass[index] = mass[index] - (DT*(mass_F[index+NY+2]-mass_F[index])/DX) - (DT*(mass_G[index+1]-mass_G[index])/DY);
                                momentum_X[index] = momentum_X[index] - (DT*(momentum_F_X[index+NY+2]-momentum_F_X[index])/DX)-(DT*(momentum_G_X[index+1]-momentum_G_X[index])/DY);
                                momentum_Y[index] = momentum_Y[index] - (DT*(momentum_F_Y[index+NY+2]-momentum_F_Y[index])/DX)-(DT*(momentum_G_Y[index+1]-momentum_G_Y[index])/DY);
//                      printf("%f,%f\n",mass[i],momentum[i]);
                        }
                }
                #pragma omp for
                for(int i = 0;i < NX+2;i++){
                        for (int j = 0; j < NY+2; j++){
                                int index = i*(NY+2)+j;
                                h[index] = mass[index];
                                u[index] = momentum_X[index]/mass[index];
                                v[index] = momentum_Y[index]/mass[index];
                        }
                }
}
