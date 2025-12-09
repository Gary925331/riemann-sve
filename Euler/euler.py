import matplotlib
matplotlib.use('TkAgg') # Or 'Qt5Agg'
import matplotlib.pyplot as plt
import array as arr

DX = 0.01
X = 1
Y = 1

X_ARRAY = arr.array('f', [1])
Y_ARRAY = arr.array('f', [1])


for i in range(50):
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
    print(f"REAL_Y = {REAL_Y}")

plt.plot(X_ARRAY, Y_ARRAY)
plt.xlabel('X-axis Label')
plt.ylabel('Y-axis Label')
plt.title('My Simple Line Graph')
plt.show()
