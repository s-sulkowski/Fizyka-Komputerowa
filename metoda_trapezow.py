import matplotlib.pyplot as plt
import numpy as np
y0 = 1
lamb = -1
tstart = 0
tend = 5
delta_t = [0.01, 0.1, 1.0]
fig1, ax = plt.subplots(1,2,figsize=(10, 6))
t_analityczne = np.linspace(tstart,tend,500)
def y(t):
    return np.exp(lamb*t)
ax[0].plot(t_analityczne,y(t_analityczne),color= "green")
for dt in delta_t:
    t = np.arange(tstart,tend+dt/2,dt)
    y_trap = np.zeros(len(t))
    y_trap[0] = y0
    for n in range(len(t)-1):
        k1 = y_trap[n]*lamb
        k2 = lamb*(y_trap[n]+dt*k1)
        y_trap[n+1] = y_trap[n] + (dt/2)*(k1+k2)
    ax[0].grid(True)
    ax[1].grid(True)
   
    ax[0].plot(t,y_trap,label = 'metoda trapezow',color = 'blue',linestyle = '--')
    ax[1].plot(t,y_trap-y(t),color = 'green',label = 'blad globalny')
ax[0].legend()
ax[1].legend()

plt.show()
