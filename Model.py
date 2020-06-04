import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def dN_dt(N,t=0):
   return np.array([-kp1*N[0]*N[2]+km1*N[3]-kp1*N[0]*N[4]+km1*N[5],
                    -kp2*N[1]*N[2]+km2*N[4]-kp2*N[1]*N[3]+km2*N[5],
                    -kp1*N[0]*N[2]+km1*N[3]-kp2*N[1]*N[2]+km2*N[4],
                     kp1*N[0]*N[2]-km1*N[3]-kp2*N[1]*N[3]+km2*N[5],
                     kp2*N[1]*N[2]-km2*N[4]-kp1*N[0]*N[4]+km1*N[5]+kp3*N[5],
                     kp2*N[1]*N[3]-km2*N[5]+kp1*N[0]*N[4]-km1*N[5]-kp3*N[5],
                     kp3*N[5]-kdeg*N[6]])

N0 = np.array((pow(10,3),pow(10,4),250,0,0,0,0))

t = np.linspace(0.0, 60, 1000)

kd1 = 500
kd2 = 5000

km1 = pow(10,-1)
km2 = pow(10,-1)

kp1 = (km1/kd1)#*60
kp2 = (km2/kd2)#*60
kp3 = pow(10,2)

kdeg = pow(10,-1)

N, infodict = integrate.odeint(dN_dt, N0, t, full_output = 1)
n1, n2, n3, n4, n5, n6, n7 = N.T

f = plt.figure(figsize=(8,4),dpi=150)
ax = f.add_subplot(111)
ax.plot(t, n1, label='POI', color='red')
#ax.plot(t, n2, label='E3')
#ax.plot(t, n3, label='PROTAC')
#ax.plot(t, n4, label='POI - PROTAC')
#ax.plot(t, n5, label='E3 - PROTAC')
#ax.plot(t, n6, label='POI - PROTAC - E3')
ax.plot(t, n7, label='ub-POI', color='blue')
ax.grid('on')

plt.legend(loc=1,title='Species')
plt.title('Timecourse of the degradation of the POI')
plt.ylabel('Concentration (nM)')
plt.xlabel('Time (mins)')
plt.show()
