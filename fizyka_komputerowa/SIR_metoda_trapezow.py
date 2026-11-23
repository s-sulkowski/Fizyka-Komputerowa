import matplotlib.pyplot as plt
import numpy as np
i0 = 0.1
s0 = 0.9
r0 = 0.0
tmax = 350
dt_tab = [0.1,10]
eps = 1e-6
MAX_ITER = 150
beta = 0.06
gamma = 0.07
for dt in dt_tab:
    t = np.arange(0,tmax,dt)
    ilosc_krokow = len(t)
    iteracje = np.zeros(ilosc_krokow)
    s = np.zeros(ilosc_krokow)
    i = np.zeros(ilosc_krokow)
    r = np.zeros(ilosc_krokow)
    s[0] = s0
    i[0] = i0
    r[0] = r0
    for n in range(0,ilosc_krokow-1): 
        s[n+1] = s[n]
        i[n+1] = i[n]
        for mi in range(MAX_ITER):
            krok_czasowy = mi+1
            s_temp = s[n+1]
            i_temp = i[n+1]
            s[n+1] = s[n] - dt/2*(beta*s[n]*i[n]+beta*s_temp*i_temp)
            i[n+1] = i[n] + dt/2*(beta*s[n]*i[n]- gamma*i[n] + beta*s_temp*i_temp - gamma*i_temp)
            if(abs(s_temp-s[n+1])<eps and abs(i_temp-i[n+1])<eps):
                break
        iteracje[n+1] = krok_czasowy
        r[n+1] = 1.0 - s[n+1] - i[n+1]
    if dt==10 :
        plt.subplot(1, 3, 1)
        plt.plot(t,s,'o',label='s',color = 'red',linestyle='none')
        plt.plot(t,i,'o',label='i',color = 'green',linestyle='none')
        plt.plot(t,r,'o',label='r',color = 'blue',linestyle='none')
        plt.grid()
        plt.legend()
    else:
        plt.subplot(1, 3, 1)
        plt.plot(t,s,label='s',color = 'red')
        plt.plot(t,i,label='i',color = 'green')
        plt.plot(t,r,label='r',color = 'blue')
        plt.grid()
        plt.legend()
    if dt==0.1:
        plt.subplot(1,3,2)
        plt.plot(t,iteracje,label='dt=0.1')
        plt.grid()
        plt.legend()
    else:
        plt.subplot(1,3,3)
        plt.plot(t,iteracje,label='dt=10')
        plt.grid()
        plt.legend()
plt.show()