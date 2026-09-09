import numpy as np
import math

NX = 1000        #  /* number of X cells */
NY = 1000         # /* number of Y cells */
N = (NX+2)*(NY+2)
NIF_X = (NX+1)     # /* number of X interfaces */
NIF_Y = (NY+1)     # /* number of Y interfaces */
NIF = (NIF_X+2)*(NIF_Y+2)
ALPHA = 1.0    #/* advection speed */
L = 200.0       #   /* domain length */
D = 200.0        #  /* domain length */
DX = (L / NX)    #/* cell size */
DY = (D / NY)    #/* cell size */
MAX_TIMESTEPS = 10
T_FINAL = 7.2
g = 9.81
CFL = 0.1

u = np.zeros(N)
v = np.zeros(N)
s = np.zeros(N)
mass_F = np.zeros(NIF)
momentum_F_X = np.zeros(NIF)
momentum_F_Y = np.zeros(NIF)
mass_G = np.zeros(NIF)
momentum_G_X = np.zeros(NIF)
momentum_G_Y = np.zeros(NIF)
mass = np.zeros(NIF)
momentum_X = np.zeros(NIF)
momentum_Y = np.zeros(NIF)
h = np.zeros(N)
h_slope_X = np.zeros(N)
u_slope_X_X = np.zeros(N)
v_slope_X_Y = np.zeros(N)
h_slope_Y = np.zeros(N)
u_slope_Y_X = np.zeros(N)
v_slope_Y_Y = np.zeros(N)

for i in range(0,NX+2):
    for j in range(0,NY+2):
        index = i*(NY+2)+j
        if(i < (NX+2)/2):
            h[index] = 10;
        else:
            h[index] = 5;
        #print(i,h[i])
        u[index] = 0
        v[index] = 0
        mass[index] = h[index]
        momentum_X[index] = h[index]*u[index]
        momentum_Y[index] = h[index]*v[index]

        h_slope_X[index] = 0
        u_slope_X_X[index] = 0
        v_slope_X_Y[index] = 0
        h_slope_Y[index] = 0
        u_slope_Y_X[index] = 0
        v_slope_Y_Y[index] = 0

time = 0
for timesteps in range (0,MAX_TIMESTEPS):
    Smax_X = 0.0;
    Smax_Y = 0.0;
    for i in range(1,NX+1):
        for j in range(1,NY+1):
            index = i*(NY+2)+j
            S_R = np.abs(u[index]) + np.sqrt(g*h[index])
            if (S_R > Smax_X):
                        Smax_X = S_R;
            S_T = np.abs(v[index]) + np.sqrt(g*mass[index])
            if (S_T > Smax_Y):
                        Smax_Y = S_T;
    term_X = Smax_X / DX
    term_Y = Smax_Y / DY
    if (term_X > term_Y):
            max_term = term_X
    else:
            max_term = term_Y
    DT = CFL / max_term
    print(DT)
    if (time > T_FINAL):
                break;
    for j in range (0,NY+2):
                mass[0*(NY+2)+j] = mass[1*(NY+2)+j]
                momentum_X[0*(NY+2)+j] = momentum_X[1*(NY+2)+j] 
                momentum_Y[0*(NY+2)+j] = momentum_Y[1*(NY+2)+j]
            
                mass[(NX+1)*(NY+2)+j] = mass[NX*(NY+2)+j]
                momentum_X[(NX+1)*(NY+2)+j] = momentum_X[NX*(NY+2)+j] 
                momentum_Y[(NX+1)*(NY+2)+j] = momentum_Y[NX*(NY+2)+j]
    for i in range(0,NX+2):
        mass[i*(NY+2)+0] = mass[i*(NY+2)+1]
        momentum_X[i*(NY+2)+0] = momentum_X[i*(NY+2)+1]
        momentum_Y[i*(NY+2)+0] = -momentum_Y[i*(NY+2)+1] #// 下牆反彈
            
        mass[i*(NY+2)+NY+1] = mass[i*(NY+2)+NY]
        momentum_X[i*(NY+2)+NY+1] = momentum_X[i*(NY+2)+NY]
        momentum_Y[i*(NY+2)+NY+1] = -momentum_Y[i*(NY+2)+NY] # // 上牆反彈

print("Hello world")

