import matplotlib.pyplot as plt
import numpy as np
y0 = 1
lamb = -1
tstart = 0
tend = 5
delta_t = [0.01, 0.1, 1.0]
fig1, ax = plt.subplots(1,2,figsize=(10, 6))
t_anlityczne = np.linspace(tstart,tend,500)
def y_dok (t):
    return np.exp(lamb*t)
ax[0].plot(t_anlityczne,y_dok(t_anlityczne),label='porownanie metod', color='green')
for dt in delta_t:
    t_rk4 = np.arange(tstart,tend + dt/2,dt)
    y_rk4 = np.zeros(len(t_rk4))
    y_rk4[0] = y0
    for i in range(len(y_rk4)-1):
        y_n = y_rk4[i]
        k1 = lamb*y_n
        k2 = lamb*(y_n+((dt/2)*k1))
        k3 = lamb*(y_n+((dt/2)*k2))
        k4 = lamb*(y_n+dt*k3)
        y_rk4[i+1] = y_n+(dt/6)*(k1+2*k2+2*k3+k4)
    ax[0].grid(True)
    ax[1].grid(True)
    ax[0].plot(t_rk4,y_rk4,label='porownanie metod',color= 'red',linestyle='--')
    ax[1].plot(t_rk4,abs(y_rk4-y_dok(t_rk4)),label= 'blad pomiarowy',color = 'green')
ax[0].legend()
ax[1].legend()
plt.show()