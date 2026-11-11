#include <stdlib.h>
#include <math.h>
#include <omp.h>

void ghost(float *u,float *v,float *a,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,
float *momentum_Y,float *h,int NX,int NY,float g){	

		#pragma omp for
                for(int j = 0; j < NY+2; j++){
                        mass[0*(NY+2)+j] = mass[1*(NY+2)+j];
                        momentum_X[0*(NY+2)+j] = momentum_X[1*(NY+2)+j]; 
                        momentum_Y[0*(NY+2)+j] = momentum_Y[1*(NY+2)+j];

                        mass[(NX+1)*(NY+2)+j] = mass[NX*(NY+2)+j];
                        momentum_X[(NX+1)*(NY+2)+j] = momentum_X[NX*(NY+2)+j]; 
                        momentum_Y[(NX+1)*(NY+2)+j] = momentum_Y[NX*(NY+2)+j];
			
			int L = 0*(NY+2)+j;
                        h[L] = mass[L];
                        u[L] = momentum_X[L]/mass[L];
                        v[L] = momentum_Y[L]/mass[L];
                        a[L] = sqrtf(g*mass[L]);

                        int R = (NX+1)*(NY+2)+j;
                        h[R] = mass[R];
                        u[R] = momentum_X[R]/mass[R];
                        v[R] = momentum_Y[R]/mass[R];
                        a[R] = sqrtf(g*mass[R]);
                }
                #pragma omp for
                for(int i = 0; i < NX+2; i++){
                        mass[i*(NY+2)+0] = mass[i*(NY+2)+1];
                        momentum_X[i*(NY+2)+0] = momentum_X[i*(NY+2)+1];
                        momentum_Y[i*(NY+2)+0] = -momentum_Y[i*(NY+2)+1]; // 下牆反彈

                        mass[i*(NY+2)+NY+1] = mass[i*(NY+2)+NY];
                        momentum_X[i*(NY+2)+NY+1] = momentum_X[i*(NY+2)+NY];
                        momentum_Y[i*(NY+2)+NY+1] = -momentum_Y[i*(NY+2)+NY]; // 上牆反彈
			
			int B = i*(NY+2)+0;
                        h[B] = mass[B];
                        u[B] = momentum_X[B]/mass[B];
                        v[B] = momentum_Y[B]/mass[B];
                        a[B] = sqrtf(g*mass[B]);

                        int T = i*(NY+2)+NY+1;
                        h[T] = mass[T];
                        u[T] = momentum_X[T]/mass[T];
                        v[T] = momentum_Y[T]/mass[T];
                        a[T] = sqrtf(g*mass[T]);
                }
/*                #pragma omp for
                for(int i = 0; i < NX+2; i++){
                        for (int j = 0; j < NY+2; j++){
                                int index = i*(NY+2)+j;
                                h[index] = mass[index];
                                u[index] = momentum_X[index]/mass[index];
                                v[index] = momentum_Y[index]/mass[index];
				a[index] = sqrtf(g*mass[index]);
                        }
                }*/
}
