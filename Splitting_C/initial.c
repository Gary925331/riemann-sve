#include <stdlib.h>
#include <math.h>

void initial(float *u,float *v,float *s,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,
float *momentum_Y,float *h,float *h_slope_X,float *u_slope_X_X,float *v_slope_X_Y,float *h_slope_Y,float *u_slope_Y_X,float *v_slope_Y_Y){
        for (int i = 0;i < NX+2;i++){
                for (int j = 0;j < NY+2;j++){
                        int index = i*(NY+2)+j;
                        if(i < (NX+2)/2){
                                h[index] = 10;
                        }else{
                                h[index] = 5;
                        }
//              printf("h[%d]=%f\n",i,h[i]);
                        u[index] = 0;
                        v[index] = 0;
                        mass[index] = h[index];
                        momentum_X[index] = h[index]*u[index];
                        momentum_Y[index] = h[index]*v[index];

                        h_slope_X[index] = 0;
                        u_slope_X_X[index] = 0;
                        v_slope_X_Y[index] = 0;
                        h_slope_Y[index] = 0;
                        u_slope_Y_X[index] = 0;
                        v_slope_Y_Y[index] = 0;
                }
        }
}
