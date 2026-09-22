#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "flux.h"

void Calculation(float *u,float *v,float *s,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,
float *momentum_Y,float *h,float *h_slope_X,float *u_slope_X_X,float *v_slope_X_Y,float *h_slope_Y,float *u_slope_Y_X,float *v_slope_Y_Y,int type,int NX,int NY,float DX,float DY,int g){	
	#pragma omp for	
	for (int i = 1; i < NX+2; i++){
			for (int j = 0; j < NY+2; j++){
				int index = i*(NY+2)+j;
            				float mass_l = mass[index-NY-2] + 0.5*DX*h_slope_X[index-NY-2];
            				float mass_r = mass[index] - 0.5*DX*h_slope_X[index];
            				float u_X_l = u[index-NY-2] + 0.5*DX*u_slope_X_X[index-NY-2];
            				float u_X_r = u[index] - 0.5*DX*u_slope_X_X[index];
					float v_Y_l = v[index-NY-2] + 0.5*DX*v_slope_X_Y[index-NY-2];
                                	float v_Y_r = v[index] - 0.5*DX*v_slope_X_Y[index];
					if (j <= NY*96/200 || j >= NY*171/200) {
            					if (i == (NX/2)+1) {
                					mass_r = mass_l;
                					u_X_r = -u_X_l;
                					v_Y_r = v_Y_l;
            					}
            					else if (i == (NX/2)+2) {
                					mass_l = mass_r;
                					u_X_l = -u_X_r;
                					v_Y_l = v_Y_r;
            					}
        				}
            				float u_l = u_X_l;
            				float u_r = u_X_r;
					float v_l = v_Y_l;
					float v_r = v_Y_r;
            				float mass_left = mass_l * u_l;
            				float mass_right = mass_r * u_r;
            				float mom_left_X = mass_l * u_l * u_l + 0.5*g*mass_l*mass_l;
            				float mom_right_X = mass_r * u_r * u_r + 0.5*g*mass_r*mass_r;
					float mom_left_Y = mass_l * v_l * u_l;
                                	float mom_right_Y = mass_r * v_r * u_r; 
					
					if (type == FLUX_HLL) {
            					float SL1 = (u_l) - sqrt(g*mass_l);
            					float SR1 = (u_r) - sqrt(g*mass_r);
                                		float S1; //SL
						if (SL1 > SR1) {
    							S1 = SR1;
						} else {
    							S1 = SL1;
						}
						float SL2 = (u_l) + sqrt(g*mass_l);
						float SR2 = (u_r) + sqrt(g*mass_r);
						float S2; //SR
						if (SL2 > SR2) {
                                                	S2 = SL2;
                                        	} else {
                                                	S2 = SR2;
                                        	}
						if(S1 > 0){
							mass_F[index] = mass_left;
							momentum_F_X[index] = mom_left_X;
							momentum_F_Y[index] = mom_left_Y;
						}else if(S2 < 0){
							mass_F[index] = mass_right;
                                                	momentum_F_X[index] = mom_right_X;
                                                	momentum_F_Y[index] = mom_right_Y;
						}else{
            						mass_F[index] = (S2*mass_left - S1*mass_right)/(S2-S1) + S1*S2*(mass_r - mass_l)/(S2-S1);
            						momentum_F_X[index] = (S2*mom_left_X - S1*mom_right_X)/(S2-S1) + S1*S2*(mass_r * u_r - mass_l * u_l)/(S2-S1);
							momentum_F_Y[index] = (S2*mom_left_Y - S1*mom_right_Y)/(S2-S1) + S1*S2*(mass_r * v_r - mass_l * v_l)/(S2-S1);
						}
					}else{
						float S_L = fabs(u_l) + sqrt(g*mass_l);
                                        	float S_R = fabs(u_r) + sqrt(g*mass_r);
                                        	float S;
                                        	if (S_L > S_R) {
                                                	S = S_L;
                                        	} else {
                                                	S = S_R;
                                        	}
                                        	mass_F[index] = 0.5*(mass_left + mass_right) - 0.5*S*(mass_r - mass_l);
                                        	momentum_F_X[index] = 0.5*(mom_left_X + mom_right_X) - 0.5*S*(mass_r * u_r - mass_l * u_l);
                                        	momentum_F_Y[index] = 0.5*(mom_left_Y + mom_right_Y) - 0.5*S*(mass_r * v_r - mass_l * v_l);
					}

			}
        	}
		#pragma omp for
		//Y direction flux
		for (int i = 0; i < NX+2; i++){
                        for (int j = 1; j < NY+2; j++){
                                int index = i*(NY+2)+j;
				//int index1 = i*(NIF_Y+2)+j;
                                float mass_B = mass[index-1] + 0.5*DY*h_slope_Y[index-1];
                                float mass_T = mass[index] - 0.5*DY*h_slope_Y[index];
                                float u_X_B = u[index-1] + 0.5*DY*u_slope_Y_X[index-1];
                                float u_X_T = u[index] - 0.5*DY*u_slope_Y_X[index];
				float v_Y_B = v[index-1] + 0.5*DY*v_slope_Y_Y[index-1];
                                float v_Y_T = v[index] - 0.5*DY*v_slope_Y_Y[index];
				if(i == (NX+2)/2){
					if(j == NY*96/200 +1){
						mass_B = mass_T;
						u_X_B = u_X_T;
						v_Y_B = -v_Y_T;
					}
					if(j == NY*171/200){
						mass_T = mass_B;
						u_X_T = u_X_B;
						v_Y_T = -v_Y_B;
					}
				}
                                float u_B = u_X_B;
                                float u_T = u_X_T;
				float v_B = v_Y_B;
				float v_T = v_Y_T;
                                float mass_Bottom = mass_B * v_B;
                                float mass_Top = mass_T * v_T;
                                float mom_Bottom_X = mass_B * u_B * v_B;
                                float mom_Top_X = mass_T * u_T * v_T;
				float mom_Bottom_Y = mass_B * v_B * v_B + 0.5*g*mass_B*mass_B;
                                float mom_Top_Y = mass_T * v_T * v_T + 0.5*g*mass_T*mass_T;
			
				if (type == FLUX_HLL) {
                                	float SB1 = (v_B) - sqrt(g*mass_B);
                                	float ST1 = (v_T) - sqrt(g*mass_T);
                                	float S1; //SL
                                	if (ST1 > SB1){
                                        	S1 = SB1;
                                	}else{
						S1 = ST1;
					}
					float SB2 = (v_B) + sqrt(g*mass_B);
                                	float ST2 = (v_T) + sqrt(g*mass_T);
                                	float S2; //SR
                                	if (ST2 > SB2){
                                        	S2 = ST2;
                                	}else{
                                        	S2 = SB2;
                                	}
					if(S1 > 0){
						mass_G[index] = mass_Bottom;
						momentum_G_X[index] = mom_Bottom_X;
						momentum_G_Y[index] = mom_Bottom_Y;
					}else if(S2 < 0){
						mass_G[index] = mass_Top;
                                        	momentum_G_X[index] = mom_Top_X;
                                        	momentum_G_Y[index] = mom_Top_Y;
					}else{
                                		mass_G[index] = (S2*mass_Bottom - S1*mass_Top)/(S2-S1) + S1*S2*(mass_T - mass_B)/(S2-S1);
                                		momentum_G_X[index] = (S2*mom_Bottom_X - S1*mom_Top_X)/(S2-S1) + S1*S2*(mass_T * u_T - mass_B * u_B)/(S2-S1);
						momentum_G_Y[index] = (S2*mom_Bottom_Y - S1*mom_Top_Y)/(S2-S1) + S1*S2*(mass_T * v_T - mass_B * v_B)/(S2-S1);
					}
				}else{
					float S_B = fabs(v_B) + sqrt(g*mass_B);
                                	float S_T = fabs(v_T) + sqrt(g*mass_T);
                                	float S;
                                	if (S_T > S_B){
                                        	S = S_T;
                                	}else{
                                        	S = S_B;
                                	}

                                	mass_G[index] = 0.5*(mass_Bottom + mass_Top) - 0.5*S*(mass_T - mass_B);
                                	momentum_G_X[index] = 0.5*(mom_Bottom_X + mom_Top_X) - 0.5*S*(mass_T * u_T - mass_B * u_B);
                                	momentum_G_Y[index] = 0.5*(mom_Bottom_Y + mom_Top_Y) - 0.5*S*(mass_T * v_T - mass_B * v_B);
				}
                        }
                }
}
