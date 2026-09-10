import numpy as np
import math
import matplotlib.pyplot as plt
from minmod import minmod
from flux import flux
from state import state

NX = 100        #  /* number of X cells */
NY = 100         # /* number of Y cells */
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

for i in range(0,NX+2):
	for j in range(0,NY+2):
		index = i*(NY+2)+j
		if(i < (NX+2)/2):
			h[index] = 10
		else:
			h[index] = 5
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
	minmod(h, u, v, NX, NY, DX, DY, h_slope_X, h_slope_Y, u_slope_X_X, u_slope_Y_X, v_slope_X_Y, v_slope_Y_Y)
	flux(mass, u, v, NX, NY, DX, DY, g, h_slope_X, h_slope_Y, u_slope_X_X, u_slope_Y_X, v_slope_X_Y, v_slope_Y_Y, mass_F, momentum_F_X, momentum_F_Y, mass_G, momentum_G_X, momentum_G_Y)
	state(NX, NY, DX, DY, DT, mass, momentum_X, momentum_Y, mass_F, momentum_F_X, momentum_F_Y, mass_G, momentum_G_X, momentum_G_Y, h, u, v) 
	time += DT
print("Hello world")
X_grid = np.zeros((NX, NY))
Y_grid = np.zeros((NX, NY))
Z_grid = np.zeros((NX, NY))  # 這裡以 mass (水深) 為例

# 2. 把 1D 陣列裡面的「內部網格」資料萃取出來，轉成 2D 矩陣
for j in range(1, NX+1):
    for k in range(1, NY+1):
        index = j * (NY + 2) + k
        
        # 存入 2D 矩陣中 (Python 矩陣索引從 0 開始，所以是 j-1, k-1)
        X_grid[j-1, k-1] = (j + 0.5) * DX
        Y_grid[j-1, k-1] = (k + 0.5) * DY
        Z_grid[j-1, k-1] = mass[index]  # 如果想畫 u 速度，可以改成 u[index]

# 3. 開始畫 3D 圖
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 使用 plot_surface 畫出曲面圖 (類似 MATLAB 的 surf / mesh)
# cmap='viridis' 會幫水深加上漂亮的漸層顏色
surf = ax.plot_surface(X_grid, Y_grid, Z_grid, cmap='viridis', edgecolor='none')

# 設定座標軸標籤
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Water Depth / mass (m)')
ax.set_title('2D Shallow Water Equations - Water Depth')

# 加上顏色條 (Colorbar)
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label='Depth')

# 顯示視窗 (這行執行後會跳出視窗，你可以用滑鼠旋轉 3D 圖！)
plt.show()
