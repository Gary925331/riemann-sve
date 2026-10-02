#ifndef MCLIMITER_H
#define MCLIMITER_H
 
#define MINMOD 0
#define MC     1

void mclimiter(float *u,float *v,float *s,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,
float *momentum_Y,float *h,float *h_slope_X,float *u_slope_X_X,float *v_slope_X_Y,float *h_slope_Y,float *u_slope_Y_X,float *v_slope_Y_Y,int NX,int NY,float DX,float DY);

#endif
