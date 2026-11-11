#ifndef FLUX_H
#define FLUX_H

#define FLUX_HLL 0
#define FLUX_RUS  1

void Calculation(float *u,float *v,float *a,float *mass_F,float *momentum_F_X,float *momentum_F_Y,float *mass_G,float *momentum_G_X,float *momentum_G_Y,float *mass,float *momentum_X,
float *momentum_Y,float *h,int type,int type_limiter,int NX,int NY,float DX,float DY,float g);

#endif
