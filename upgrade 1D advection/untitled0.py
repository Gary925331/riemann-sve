import numpy as np
import matplotlib.pyplot as plt

# 讀取 fluxes.dat 檔案
data = np.loadtxt("fluxes.dat")

# 分成兩個欄位：index 與 flux
j = data[:, 0]
F = data[:, 1]

# 繪圖
plt.figure(figsize=(8, 4))
plt.plot(j, F, 'b-', linewidth=2, label='Flux $F_j = \\alpha u_{interface}$')
plt.xlabel('Interface index j')
plt.ylabel('Flux F')
plt.title('Computed Flux Distribution (α = 0.1)')
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()
