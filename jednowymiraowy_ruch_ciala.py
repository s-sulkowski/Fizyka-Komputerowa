import matplotlib.pyplot as plt
import numpy as np
k = 1
m = 1
x0 = 0
v0 = 0
dt = 0.0001
omega = 1
Omega = 1
gamma = 0
f0 = 0.5
t0 = 2*(np.pi)/omega
t = np.arange(0,4*t0,dt)
beta = 0
x = np.zeros(len(t))
v = np.zeros(len(t))
fig1, ax = plt.subplots(1,2,figsize=(18, 10))
x[0] = x0
v[0] = v0
def F(t):
    return f0*np.cos(Omega*t+gamma)
for n in range(len(t)-1): 
    k1x = v[n]
    k1v = -pow(omega,2)*x[n] - beta*v[n] + F(t[n])
    k2x = v[n] + (dt/2)*k1v
    k2v = -pow(omega,2)*(x[n]+(dt/2)*k1x) - beta*(v[n]+(dt/2)*k1v) + F(t[n]+(1/2)*dt)
    k3x = v[n] + (dt/2)*k2v
    k3v = -pow(omega,2)*(x[n]+(dt/2)*k2x) - beta*(v[n]+(dt/2)*k2v) + F(t[n]+(1/2)*dt)
    k4x =  v[n] + dt*k3v
    k4v = -pow(omega,2)*(x[n]+dt*k3x) - beta*(v[n]+dt*k3v) + F(t[n]+dt)
    x[n+1] = x[n] + (dt/6)*(k1x+2*k2x + 2*k3x+k4x)
    v[n+1] = v[n] + (dt/6)*(k1v+2*k2v + 2*k3v+k4v)
ax[0].grid(True)
ax[1].grid(True)
ax[0].plot(t,x,label = "x(t)",color = "green")
ax[0].plot(t,v,label = "v(t)",color = "red")
ax[0].legend()
ax[1].plot(t,0.5*k*pow(x,2),label = "Epot",color="orange")
ax[1].plot(t,0.5*m*pow(v,2),label = "Ekin",color = "green")
ax[1].plot(t,0.5*k*pow(x,2)+0.5*m*pow(v,2),label = "Ecalk",color = "red")
ax[1].legend()
plt.show()