#include <stdlib.h>
#include <math.h>
#include <omp.h>
void state(float *restrict u,float *restrict v,float *restrict a,float *restrict mass_F,float *restrict momentum_F_X,float *restrict momentum_F_Y,float *restrict mass_G,float *restrict momentum_G_X,float *restrict momentum_G_Y,float *restrict mass,float *restrict momentum_X,
float *restrict momentum_Y,float *restrict h,int NX,int NY,float DX,float DY,float DT,float g,float *Smax_X,float *Smax_Y){
		float sx = 0.0f;
		float sy = 0.0f;
		#pragma omp for
                for(int i = 1;i < NX+1;i++){
                        for (int j = 1; j < NY+1; j++){
                                int index = i*(NY+2)+j;
                                //int index1 = i*(NIF_Y+2)+j;
                                if (i == NX/2 +1 && (j <= NY*96/200 || j >= NY*171/200)) {
                                        float S_R = fabsf(u[index]) + a[index];
                                        if (S_R > sx) sx = S_R;
                                        float S_T = fabsf(v[index]) + a[index];
                                        if (S_T > sy) sy = S_T;
					continue; 
                                }
                                mass[index] = mass[index] - (DT*(mass_F[index+NY+2]-mass_F[index])/DX) - (DT*(mass_G[index+1]-mass_G[index])/DY);
                                momentum_X[index] = momentum_X[index] - (DT*(momentum_F_X[index+NY+2]-momentum_F_X[index])/DX)-(DT*(momentum_G_X[index+1]-momentum_G_X[index])/DY);
                                momentum_Y[index] = momentum_Y[index] - (DT*(momentum_F_Y[index+NY+2]-momentum_F_Y[index])/DX)-(DT*(momentum_G_Y[index+1]-momentum_G_Y[index])/DY);
//                      printf("%f,%f\n",mass[i],momentum[i]);
				float inv = 1.0f/mass[index];
				h[index] = mass[index];
                                u[index] = momentum_X[index]*inv;
                                v[index] = momentum_Y[index]*inv;
                                a[index] = sqrtf(g*mass[index]);
				float S_R = fabsf(u[index]) + a[index];
                                if (S_R > sx) sx = S_R;
                                float S_T = fabsf(v[index]) + a[index];
                                if (S_T > sy) sy = S_T; 
                        }
                }
		*Smax_X = sx;      
                *Smax_Y = sy;
               /* #pragma omp for
                for(int i = 0;i < NX+2;i++){
                        for (int j = 0; j < NY+2; j++){
                                int index = i*(NY+2)+j;
				float inv = 1.0f/mass[index];
                                h[index] = mass[index];
                                u[index] = momentum_X[index]*inv;
                                v[index] = momentum_Y[index]*inv;
				a[index] = sqrtf(g*mass[index]); 
                        }
                }*/
}
