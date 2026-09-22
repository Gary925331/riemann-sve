#ifndef FLUX_H
#define FLUX_H

#define FLUX_HLL 0
#define FLUX_RUS  1

void Calculation(float *u,float *v,float *s,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,
float *momentum_Y,float *h,float *h_slope_X,float *u_slope_X_X,float *v_slope_X_Y,float *h_slope_Y,float *u_slope_Y_X,float *v_slope_Y_Y,int type,int NX,int NY,float DX,float DY,int g);

#endif
