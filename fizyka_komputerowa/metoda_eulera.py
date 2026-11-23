import matplotlib.pyplot as plt
import numpy as np
fig1, ax  = plt.subplots(1,2,figsize=(10,6))
y0 = 1
lamb = -1
delta_t = [0.01,0.1,1]
t0 = 0
tend = 5
x_dok = np.linspace(t0,tend,500)
def y_dok (t):
    return np.exp(lamb * t )

ax[0].plot(x_dok,y_dok(x_dok),label='analityczna',color='red')
for dt in delta_t:
    t_eu = np.arange(t0,tend+dt/2,dt)
    y_eu = np.zeros(len(t_eu))
    y_eu[0] = y0
    for i in range(len(y_eu)-1):
        y_n = y_eu[i] 
        y_eu[i+1] = y_n + (dt* lamb * y_n)
    ax[0].grid(True)
    ax[1].grid(True)
    ax[0].legend()
    ax[0].plot(t_eu,y_eu,label = 'eulera',color = 'orange',linestyle="--")
    ax[1].legend()
    ax[1].plot(t_eu,y_eu-y_dok(t_eu),label = 'blad',color = 'green')
plt.show()