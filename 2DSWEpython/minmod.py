import numpy as np

def minmod(h,u,v):
	for i in range (1,NX+1):
		for j in range (0,NY+2):
			index = i*(NY+2)+j
			forward = (h[index+NY+2] - h[index])/DX
			backward = (h[index] - h[index-NY-2])/DX
			if(forward*backward < 0):
				h_slope_X[index] = 0
			else:
				if(np.abs(forward)<np.abs(backward)):
					h_slope_X[index] = forward
				else:
					h_slope_X[index] = backward
	for i in range (0,NX+2):
		for j in range (1,NY+1):
			index = i*(NY+2)+j
			forward = (h[index+1] - h[index])/DY
			backward = (h[index] - h[index-1])/DY
			if(forward*backward < 0):
				h_slope_Y[index] = 0
			else:
				if(fabs(forward)<fabs(backward)):
					h_slope_Y[index] = forward
				else:
					h_slope_Y[index] = backward
	for i in range (1,NX+1):
		for j in range (0,NY+2):
			index = i*(NY+2)+j
			forward = (u[index+NY+2] - u[index])/DX
			backward = (u[index] - u[index-NY-2])/DX
			if(forward*backward < 0):
				u_slope_X_X[index] = 0
			else:
				if(fabs(forward)<fabs(backward)):
					u_slope_X_X[index] = forward
				else:
					u_slope_X_X[index] = backward
	for i in range (1,NX+1):
		for j in range (0,NY+2):
			index = i*(NY+2)+j
			forward = (v[index+NY+2] - v[index])/DX
			backward = (v[index] - v[index-NY-2])/DX
			if(forward*backward < 0):
				v_slope_X_Y[index] = 0
			else:
				if(fabs(forward)<fabs(backward)):
					v_slope_X_Y[index] = forward
				else:
					v_slope_X_Y[index] = backward
	for i in range (0,NX+2):
		for j in range (1,NY+1):
			index = i*(NY+2)+j
			forward = (u[index+1] - u[index])/DY
			backward = (u[index] - u[index-1])/DY
			if(forward*backward < 0):
				u_slope_Y_X[index] = 0
			else:
				if(fabs(forward)<fabs(backward)):
					u_slope_Y_X[index] = forward
				else:
					u_slope_Y_X[index] = backward
	for i in range (0,NX+2):
		for j in range (1,NY+1):
			index = i*(NY+2)+j
			forward = (v[index+1] - v[index])/DY
			backward = (v[index] - v[index-1])/DY
			if(forward*backward < 0):
				v_slope_Y_Y[index] = 0
			else:
				if(fabs(forward)<fabs(backward)):
					v_slope_Y_Y[index] = forward
				else:
					v_slope_Y_Y[index] = backward

