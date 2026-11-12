#ifndef FLUX_H
#define FLUX_H

#define FLUX_HLL 0
#define FLUX_RUS  1

void Calculation(float *restrict u,float *restrict v,float *restrict a,float *restrict mass_F,float *restrict momentum_F_X,float *restrict momentum_F_Y,float *restrict mass_G,float *restrict momentum_G_X,float *restrict momentum_G_Y,float *restrict mass,float *restrict momentum_X,
float *restrict momentum_Y,float *restrict h,int type,int type_limiter,int NX,int NY,float DX,float DY,float g);

#endif
