import numpy as np

def flux(mass,u,v):
	#X direction
	for i in range(1,NX+2):
		for j in range (0,NY+2):
			index = i*(NY+2)+j
			mass_l = mass[index-NY-2] + 0.5*DX*h_slope_X[index-NY-2]
			mass_r = mass[index] - 0.5*DX*h_slope_X[index]
			u_X_l = u[index-NY-2] + 0.5*DX*u_slope_X_X[index-NY-2]
			u_X_r = u[index] - 0.5*DX*u_slope_X_X[index]
			v_Y_l = v[index-NY-2] + 0.5*DX*v_slope_X_Y[index-NY-2]
			v_Y_r = v[index] - 0.5*DX*v_slope_X_Y[index]
			if (j <= 96 or j >= 171):
				if (i == 101):
					mass_r = mass_l
					u_X_r = -u_X_l
					v_Y_r = v_Y_l
				elif (i == 102):
					mass_l = mass_r
					u_X_l = -u_X_r
					v_Y_l = v_Y_r

			u_l = u_X_l
			u_r = u_X_r
			v_l = v_Y_l
			v_r = v_Y_r
			mass_left = mass_l * u_l
			mass_right = mass_r * u_r
			mom_left_X = mass_l * u_l * u_l + 0.5*g*mass_l*mass_l
			mom_right_X = mass_r * u_r * u_r + 0.5*g*mass_r*mass_r
			mom_left_Y = mass_l * v_l * u_l
			mom_right_Y = mass_r * v_r * u_r 

			SL1 = (u_l) - np.sqrt(g*mass_l)
			SR1 = (u_r) - np.sqrt(g*mass_r)
			if (SL1 > SR1):
				S1 = SR1 #SL
			else:
				S1 = SL1
			SL2 = (u_l) + np.sqrt(g*mass_l)
			SR2 = (u_r) + sqrt(g*mass_r)
			if (SL2 > SR2):
				S2 = SL2 #SR
			else:
				S2 = SR2
			if(S1 > 0):
				mass_F[index] = mass_left
				momentum_F_X[index] = mom_left_X
				momentum_F_Y[index] = mom_left_Y
			elif(S2 < 0):
				mass_F[index] = mass_right
				momentum_F_X[index] = mom_right_X
				momentum_F_Y[index] = mom_right_Y
			else:
				mass_F[index] = (S2*mass_left - S1*mass_right)/(S2-S1) + S1*S2*(mass_r - mass_l)/(S2-S1)
				momentum_F_X[index] = (S2*mom_left_X - S1*mom_right_X)/(S2-S1) + S1*S2*(mass_r * u_r - mass_l * u_l)/(S2-S1)
				momentum_F_Y[index] = (S2*mom_left_Y - S1*mom_right_Y)/(S2-S1) + S1*S2*(mass_r * v_r - mass_l * v_l)/(S2-S1)

	#Y direction flux
	for i in range (0,NX+2):
		for j in range (1,NY+2):
			index = i*(NY+2)+j
			mass_B = mass[index-1] + 0.5*DY*h_slope_Y[index-1]
			mass_T = mass[index] - 0.5*DY*h_slope_Y[index]
			u_X_B = u[index-1] + 0.5*DY*u_slope_Y_X[index-1]
			u_X_T = u[index] - 0.5*DY*u_slope_Y_X[index]
			v_Y_B = v[index-1] + 0.5*DY*v_slope_Y_Y[index-1]
			v_Y_T = v[index] - 0.5*DY*v_slope_Y_Y[index]
			if(i == 101):
				if(j == 97 ):
					mass_B = mass_T
					u_X_B = u_X_T
					v_Y_B = -v_Y_T
				if(j == 171):
					mass_T = mass_B
					u_X_T = u_X_B
					v_Y_T = -v_Y_B
			u_B = u_X_B
			u_T = u_X_T
			v_B = v_Y_B
			v_T = v_Y_T
			mass_Bottom = mass_B * v_B
			mass_Top = mass_T * v_T
			mom_Bottom_X = mass_B * u_B * v_B
			mom_Top_X = mass_T * u_T * v_T
			mom_Bottom_Y = mass_B * v_B * v_B + 0.5*g*mass_B*mass_B
			mom_Top_Y = mass_T * v_T * v_T + 0.5*g*mass_T*mass_T;

			SB1 = (v_B) - sqrt(g*mass_B)
			ST1 = (v_T) - sqrt(g*mass_T)
			if (ST1 > SB1):
				S1 = SB1 #SL
			else:
				S1 = ST1
			SB2 = (v_B) + sqrt(g*mass_B)
			ST2 = (v_T) + sqrt(g*mass_T)
			if (ST2 > SB2):
				S2 = ST2
			else:
				S2 = SB2
			if(S1 > 0):
				mass_G[index] = mass_Bottom
				momentum_G_X[index] = mom_Bottom_X
				momentum_G_Y[index] = mom_Bottom_Y
			elif(S2 < 0):
				mass_G[index] = mass_Top
				momentum_G_X[index] = mom_Top_X
				momentum_G_Y[index] = mom_Top_Y
			else:
				mass_G[index] = (S2*mass_Bottom - S1*mass_Top)/(S2-S1) + S1*S2*(mass_T - mass_B)/(S2-S1)
				momentum_G_X[index] = (S2*mom_Bottom_X - S1*mom_Top_X)/(S2-S1) + S1*S2*(mass_T * u_T - mass_B * u_B)/(S2-S1)
				momentum_G_Y[index] = (S2*mom_Bottom_Y - S1*mom_Top_Y)/(S2-S1) + S1*S2*(mass_T * v_T - mass_B * v_B)/(S2-S1)
				
                        
                
