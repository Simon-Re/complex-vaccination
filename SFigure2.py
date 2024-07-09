from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

#computing the probability of survival (exitinction) using a branching process
def f(s,R0):
    return np.exp(-R0*(1-s)) - s

#evaluate the branching process for different values of R0 using scipy.fsolve
Rspace = np.linspace(0.1,100.,100000)
solspace = []
for R0 in Rspace:
    a=fsolve(f,0,args = R0)[0]
    solspace.append(1-a)
 
#Plotting
#########
    
fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(13, 7),dpi=300)
fig.tight_layout(pad=8.0)

ax1.plot(Rspace,solspace,label='Poisson Branching Process',linewidth=5)
ax1.plot(Rspace,np.heaviside(1.-1./Rspace,0.)*(1.-1./(Rspace)),linewidth=5,label='Birth Death Process') #compare with birth death process prediction
ax1.set_xlim(0,10)
ax1.tick_params(labelsize=15)
ax1.set_ylabel(r'$1 -p_{ext}$',size=20)
ax1.set_xlabel(r'$R_0$',size=20)
ax1.legend(prop={'size':15})
ax1.set_title('(a)',size=20,loc='left')

ax2.plot(solspace,np.heaviside(1.-1./Rspace,0.)*(1.-1./(Rspace)), 'r-',linewidth=5)
ax2.plot(solspace,solspace, 'k--',alpha=0.5,linewidth=5)
ax2.set_ylabel(r'$1 -p_{ext}$' +'\n (Birth Death Process)',size=20)
ax2.set_xlabel(r'$1 -p_{ext}$' +'\n (Poisson Branching Process)',size=20)
ax2.set_title('(b)',size=20,loc='left')
ax2.tick_params(labelsize=15)
plt.savefig('SFigure2.jpg')
plt.show()
