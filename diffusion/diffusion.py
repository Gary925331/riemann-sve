import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg') # Or 'Qt5Agg'
import matplotlib.pyplot as plt

L = 1
NO_CELLS = 200
U = 0.5
ALPHA = 1e-4
DX = L / NO_CELLS
T_END = 0.4
DIFFUSION_CFL = 0.001  # ALPHA*DT/(DX^2)
DT = (DIFFUSION_CFL*DX*DX)/ALPHA
NO_TIME_STEPS = int(T_END/DT)
NO_INTERFACES = NO_CELLS + 1

x = np.linspace(0.0, 1.0 ,NO_CELLS)
H_INIT = np.zeros(NO_CELLS)
H = np.zeros(NO_CELLS)

H_new = np.zeros(NO_CELLS)
F = np.zeros(NO_INTERFACES)


start = int(0.2 / DX)
end = int(0.4 / DX)
for i in range(start, end):
    H[i] = 1.0

print(f"Initial pulse in cells {start} to {end - 1}")


H_INIT[:] = H[:]

print(f"H = {H}")

for step in range(NO_TIME_STEPS):
    print(f"Computing steps {step}")

    for i in range(1, NO_INTERFACES - 1):
        #print("Computing flux at interface", i)

        F[i] = 0

        if (U > 0):
            F[i] += U * H[i - 1]
        else:
            F[i] += U * H[i]
        
        F[i] += -ALPHA*(H[i] - H[i - 1]) / DX

     

    
    for i in range(NO_CELLS):
        H_new[i] = H[i] - (DT / DX) * (F[i+1] - F[i])

    H[:] = H_new[:]

    


plt.plot(x, H)
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.savefig("diffusion_graph.png", dpi=1000, bbox_inches='tight')


