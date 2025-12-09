import matplotlib
matplotlib.use('TkAgg') # Or 'Qt5Agg'
import matplotlib.pyplot as plt
import array as arr

DX = 0.001
NO_STEPS = 500
X = 1
Y = 1

X_ARRAY = arr.array('f', [1])
Y_ARRAY = arr.array('f', [1])
REAL_Y_ARRAY = arr.array('f', [1])

for i in range(NO_STEPS):
    print(f"Computing step {i}")
    DY_DX = (Y*Y +2*X*Y)/(X*X)
    print(f"X,Y = ({X},{Y}), DY_DX = {DY_DX}")
    DY = DX*DY_DX
    Y = Y + DY
    X = X + DX
    X_ARRAY.append(X)
    Y_ARRAY.append(Y)
    print(f"(X,Y) = {X}, {Y}")
    REAL_Y = (1/2)*X*X/(1-0.5*X)
    REAL_Y_ARRAY.append(REAL_Y)
    print(f"REAL_Y = {REAL_Y}")

plt.plot(X_ARRAY, Y_ARRAY, label = 'Euler_Y')
plt.plot(X_ARRAY,REAL_Y_ARRAY, label = 'REAL_Y')
plt.xlabel('X')
plt.ylabel('Y')
title_string = f"Euler solution to dy/dx=(y*y+2*x*y)/x*x, dx={DX}, {NO_STEPS}"
plt.title(title_string)
plt.legend()
plt.savefig("Euler_method_graph.png", dpi=1000, bbox_inches='tight')
