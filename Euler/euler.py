DX = 0.1
X = 1
Y = 1

for i in range(5):
    print(f"Computing step {i}")
    DY_DX = (Y*Y +2*X*Y)/X*X
    print(f"DY_DX = {DY_DX}")
    DY = DX*DY_DX
    Y = Y + DY
    X = X + DX
    print(f"(X,Y) = {X}, {Y}")
