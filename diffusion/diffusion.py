import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg') # Or 'Qt5Agg'
import matplotlib.pyplot as plt

L = 0.1
NO_CELLS = 200
U = 0.5
ALPHA = 0.1
DX = L / NO_CELLS
T_END = 0.5
DT = DX / ALPHA
NO_TIME_STEPS = 500
NO_INTERFACES = NO_CELLS + 1

H_INIT = np.zeros(NO_CELLS)
H = np.zeros(NO_CELLS)
H_new = np.zeros(NO_CELLS)
F = np.zeros(NO_INTERFACES)


H_INIT[:] = H[:]

for step in range(NO_TIME_STEPS):
    print(f"Computing steps {step}")

    for i in range(1, NO_INTERFACES - 1):

        F[i] = 0

        if (U > 0):
            F[i] += U * H[i - 1]
        else:
            F[i] += U * H[i]
        
        F[i] += -ALPHA*(H[i] - H[i - 1]) / DX

    
    for i in range(NO_CELLS):
        H_new[i] = H[i] - (DT / DX) * (F[i+1] - F[i])

    H[:] = H_new[:]

x = np.linspace(0.0, 1.0 ,NO_CELLS)
plt.plot(x, H)
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.show()
plt.savefig("Euler_method_graph.png", dpi=1000, bbox_inches='tight')


