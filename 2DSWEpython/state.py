import numpy as np

def state(NX, NY, DX, DY, DT, mass, momentum_X, momentum_Y, mass_F, momentum_F_X, momentum_F_Y, mass_G, momentum_G_X, momentum_G_Y, h, u, v):
	for i in range (1,NX+1):
		for j in range (1,NY+1):
			index = i*(NY+2)+j
			if (i == 101 and (j <= 96 or j >= 171)):
            			continue
			mass[index] = mass[index] - (DT*(mass_F[index+NY+2]-mass_F[index])/DX) - (DT*(mass_G[index+1]-mass_G[index])/DY)
			momentum_X[index] = momentum_X[index] - (DT*(momentum_F_X[index+NY+2]-momentum_F_X[index])/DX)-(DT*(momentum_G_X[index+1]-momentum_G_X[index])/DY)
			momentum_Y[index] = momentum_Y[index] - (DT*(momentum_F_Y[index+NY+2]-momentum_F_Y[index])/DX)-(DT*(momentum_G_Y[index+1]-momentum_G_Y[index])/DY)
	for i in range (0,NX+2):
		for j in range (0,NY+2):
			index = i*(NY+2)+j
			h[index] = mass[index]
			u[index] = momentum_X[index]/mass[index]
			v[index] = momentum_Y[index]/mass[index]
