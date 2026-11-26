import matplotlib.pyplot as plt
import numpy as np
Nx = 50
Ny = 100 
delta_x= 0.1 
Vc = 10
epsilon = 1e-8

Vs = np.zeros((Nx+1,Ny+1))
x = np.arange(Nx+1) * delta_x
y = np.arange(Ny+1) * delta_x
X, Y = np.meshgrid(x, y, indexing='ij')
srodek_x = delta_x*Nx/2
srodek_y = delta_x*Ny/2
Rc = 1
maska = (X - srodek_x)**2 + (Y - srodek_y)**2 <= Rc**2
Vs[maska] = Vc
print(maska.sum())
plt.imshow(Vs.T, origin='lower', extent=[0, x[-1], 0, y[-1]])
plt.colorbar()
plt.show()

