import matplotlib.pyplot as plt
import numpy as np
N = 50
dt_t = 0.02
dt = 0.1
alpha = 1
mi = 1
lamb = -1

def f(t, y):
    return lamb*y

def RK4(f,t,y,dt):
    k1 = f(t,y)
    k2 = f(t + dt/2, y + dt/2 * k1)
    k3 = f(t + dt/2, y + dt/2 * k2)
    k4 = f(t + dt, y + dt * k3)
    return y+(dt/6)*(k1+2*k2+2*k3+k4)

