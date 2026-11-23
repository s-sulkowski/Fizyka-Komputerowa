import matplotlib.pyplot as plt
import numpy as np
import math
theta0= math.radians(179) 
S = 0.75
dt0 = 0.01
tmax = 10
p = 2
grav = 9.81
lamb = 1 
alpha0 = 0.5
alpha1 = 1000
R = 1
t = 0 

def f (t,theta,v0):
    return v0
def g (t,theta,v0):
    alpha = alpha0 + alpha1*np.exp(-lamb*t)
    return -grav/R*np.sin(theta)-alpha*v0

def metoda_trapezow(theta,v,dt,t):
    theta_new = theta
    v_new = v
    f_n = f(t, theta, v)
    g_n = g(t, theta, v)
    alpha_next = alpha0 + alpha1 * np.exp(-lamb * (t + dt))
    while True:
        F = theta_new - theta - (dt / 2.0) * (f_n + f(t + dt, theta_new, v_new))
        G = v_new - v - (dt / 2.0) * (g_n + g(t + dt, theta_new, v_new))
        a11 = 1.0
        a12 = -dt / 2.0
        a21 = (dt / 2.0) * (grav / R) * np.cos(theta_new)
        a22 = 1.0 + (dt / 2.0) * alpha_next

        A = np.array([[a11, a12], [a21, a22]])
        b = np.array([-F, -G])
        dx, dv = np.linalg.solve(A, b)

        theta_new = theta_new + dx
        v_new = v_new + dv

        if abs(dx) < 1e-10 and abs(dv) < 1e-10:
            break
    return theta_new,v_new

def Trapezy_step(theta,v,dt,t):
    x_n_1,v_n_1 = metoda_trapezow(theta,v,dt,t)
    x_n_2,v_n_2 = metoda_trapezow(x_n_1,v_n_1,dt,t+dt)
    x_long,v_long = metoda_trapezow(theta,v,2*dt,t)

    Ex = (x_n_2 - x_long) / (2**p - 1)
    Ev = (v_n_2 - v_long) / (2**p - 1)
    err = max(abs(Ex), abs(Ev))
    return x_n_2, v_n_2, err

tol = [1e-5,1e-2]
fig, ax = plt.subplots(2, 2, figsize=(12, 8))
for i in range(len(tol)):
    theta, v, t = theta0, 0.0, 0.0
    dt = dt0
    ts, thetas, vs, dts = [t], [theta], [v], [dt]
    while t < tmax:
        if t + 2*dt > tmax:
            dt = (tmax - t) / 2
        theta_new, v_new, err = Trapezy_step(theta, v, dt, t)
        err = max(err, 1e-16)
        if err <= tol[i]:
            theta, v = theta_new, v_new
            t += 2*dt
            ts.append(t); thetas.append(theta); vs.append(v); dts.append(dt)
        dt = S * dt * (tol[i]/err)**(1/(p+1))
    etykieta = f"tol = {tol[i]:.0e}"
    ax[0, 0].plot(ts, thetas, ".-", markersize=3, label=etykieta)
    ax[0, 1].plot(ts, vs, ".-", markersize=3, label=etykieta)
    ax[1, 0].semilogy(ts, dts, ".-", markersize=3, label=etykieta)
    ax[1, 1].plot(thetas, vs, ".-", markersize=3, label=etykieta)
    
ax[0, 0].set_xlabel("t [s]")
ax[0, 0].set_ylabel("θ [rad]")
ax[0, 0].set_title("Kąt w czasie")
ax[0, 0].grid(True)
ax[0, 0].legend()

ax[0, 1].set_xlabel("t [s]")
ax[0, 1].set_ylabel("v = dθ/dt [rad/s]")
ax[0, 1].set_title("Prędkość kątowa w czasie")
ax[0, 1].grid(True)
ax[0, 1].legend()

ax[1, 0].set_xlabel("t [s]")
ax[1, 0].set_ylabel("dt [s]")
ax[1, 0].set_title("Krok czasowy (adaptacyjny)")
ax[1, 0].grid(True)
ax[1, 0].legend()

ax[1, 1].set_xlabel("θ [rad]")
ax[1, 1].set_ylabel("v [rad/s]")
ax[1, 1].set_title("Portret fazowy")
ax[1, 1].grid(True)
ax[1, 1].legend()

plt.tight_layout()
plt.show()
