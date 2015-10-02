import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

def drange2(start, stop, step):
    numelements = int((stop-start)/float(step))
    for i in range(numelements+1):
            yield start + i*step

def X_t(t):
	XofT = 1
	if (t < 0 or t > 4):
		XofT = 0

	return XofT

def H_t(t):
	HofT = 2
	if (t < 0 or t > 2):
		HofT = 0

	return HofT


t_beg = -2
d_t = 0.1
t_end = 8
t = np.linspace(t_beg,t_end,800)

y = np.zeros((len(t)))


tau_beg = -10
d_tau = 0.01
tau_end = 10
       
for m in range(1,len(t)):
	
	for tau in drange2(tau_beg,tau_end,d_tau):
		
		y[m] = y[m]+X_t(tau)*H_t(m-tau)*d_tau;


plt.plot(y)
plt.axis([-2,8,0,5])
plt.ylabel('y(t)')
plt.show()

