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
MAX_TIMESTEPS = 50000
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

print("Hello world")

