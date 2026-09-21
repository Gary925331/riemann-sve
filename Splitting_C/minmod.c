#include <stdlib.h>
#include <math.h>
#include <omp.h>

void minmod(float *u,float *v,float *s,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,
float *momentum_Y,float *h,float *h_slope_X,float *u_slope_X_X,float *v_slope_X_Y,float *h_slope_Y,float *u_slope_Y_X,float *v_slope_Y_Y){		

		#pragma omp for
		for(int i = 1;i < NX+1;i++){
			for(int j = 0;j < NY+2;j++){
				int index = i*(NY+2)+j;
				float forward = (mass[index+NY+2] - mass[index])/DX;
				float backward = (mass[index] - mass[index-NY-2])/DX;
				if(forward*backward < 0){
					h_slope_X[index] = 0;
				}else{
					if(fabs(forward)<fabs(backward)){
						h_slope_X[index] = forward;
					}else{
						h_slope_X[index] = backward;
					}
				}
//			printf("mass_slope[%d] = %f\n",i,mass_slope[i]);
			}
		}
		#pragma omp for
		for(int i = 0;i < NX+2;i++){
                        for(int j = 1;j < NY+1;j++){
                                int index = i*(NY+2)+j;
                                float forward = (mass[index+1] - mass[index])/DY;
                                float backward = (mass[index] - mass[index-1])/DY;
                                if(forward*backward < 0){
                                        h_slope_Y[index] = 0;
                                }else{
                                        if(fabs(forward)<fabs(backward)){
                                                h_slope_Y[index] = forward;
                                        }else{
                                                h_slope_Y[index] = backward;
                                        }
                                }
//                      printf("mass_slope[%d] = %f\n",i,mass_slope[i]);
                        }
                }
		#pragma omp for
//		printf("mass_slope[%d] = %f\n",i,mass_slope[i]);
		for(int i = 1;i < NX+1;i++){
			for(int j = 0;j < NY+2;j++){
				int index = i*(NY+2)+j;
                        	float forward = (u[index+NY+2] - u[index])/DX;
                        	float backward = (u[index] - u[index-NY-2])/DX;
                        	if(forward*backward < 0){
                                	u_slope_X_X[index] = 0;
                        	}else{
                                	if(fabs(forward)<fabs(backward)){
                                        	u_slope_X_X[index] = forward;
                                	}else{
                                        	u_slope_X_X[index] = backward;
                                	}
				}
			}
			//printf("momentum_slope[%d] = %f\n",i,momentum_slope[i]);
                }
		#pragma omp for
		for(int i = 1;i < NX+1;i++){
			for(int j = 0;j < NY+2;j++){
                                int index = i*(NY+2)+j;
                        	float forward = (v[index+NY+2] - v[index])/DX;
                        	float backward = (v[index] - v[index-NY-2])/DX;
                        	if(forward*backward < 0){
                                	v_slope_X_Y[index] = 0;
                        	}else{
                                	if(fabs(forward)<fabs(backward)){
                                        	v_slope_X_Y[index] = forward;
                                	}else{
                                        	v_slope_X_Y[index] = backward;
                                	}
                        	}
			}
                        //printf("momentum_slope[%d] = %f\n",i,momentum_slope[i]);
                }
		#pragma omp for
		for(int i = 0;i < NX+2;i++){
                        for(int j = 1;j < NY+1;j++){
                                int index = i*(NY+2)+j;
                                float forward = (u[index+1] - u[index])/DY;
                                float backward = (u[index] - u[index-1])/DY;
                                if(forward*backward < 0){
                                        u_slope_Y_X[index] = 0;
                                }else{
                                        if(fabs(forward)<fabs(backward)){
                                                u_slope_Y_X[index] = forward;
                                        }else{
                                                u_slope_Y_X[index] = backward;
                                        }
                        	}
			}
                        //printf("momentum_slope[%d] = %f\n",i,momentum_slope[i]);
                }
		#pragma omp for
		for(int i = 1;i < NX+1;i++){
                        for(int j = 0;j < NY+2;j++){
                                int index = i*(NY+2)+j;
                                float forward = (v[index+1] - v[index])/DY;
                                float backward = (v[index] - v[index-1])/DY;
                                if(forward*backward < 0){
                                        v_slope_Y_Y[index] = 0;
                                }else{
                                        if(fabs(forward)<fabs(backward)){
                                                v_slope_Y_Y[index] = forward;
                                        }else{
                                                v_slope_Y_Y[index] = backward;
                                        }
                        	}
			}
                        //printf("momentum_slope[%d] = %f\n",i,momentum_slope[i]);
                }
}
