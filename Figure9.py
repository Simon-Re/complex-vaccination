import numpy as np
import matplotlib.pyplot as plt
from pyrasaic_tools import *
from scipy.optimize import fsolve

#setting model parameters
alphaspace = np.linspace(0,1,300)
n =2
wibes = 0.1
V = 0.99
dur = 5
ol = 0
pbt = 0
q = 10**-10
p = 10**-3
c = 0
rho = 0.0001
N = 10**6
base = 2.
rhospace = np.linspace(0.001,1.,100)
xlabel = r'$\rho$'

#initialization
minalphalist = []
mosaiclist = []
pyramidlist = []
bestlist = []
stdbestlist = []
stdmosaiclist = []
stdpyramidlist = []

#sample through different values of rho
for rho in rhospace:
    alpha_result_list = []
    std_result_list = []
    #check different values of alpha, to find the best one
    for alpha_strat in alphaspace:
        #calculating the probabilities of survival for each variant
        res = fix_all_strats(wibes, n, base, ol, 0, c, V, dur, p, q, rho, alpha_strat,N)
        fixonly = res[5]
        mut = res[3]
        em = res[1]
        
        #calculation of the pandemic burden for each variant
        R_IC = N*rho
        R_wt = (1- rho)*(alpha_strat*R0_func(base,1,0,n,pbt,c,V) + (1- 
                               alpha_strat)*R0_func(base,2,0,n,pbt,c,V)) + rho*base
                
        #derivation of the pandemic size of each variant
        pandemic_size = [1.-fsolve(f,0,args = R_wt)[0]]
        for i in range(1,n+1):
            R = (1- rho)*(alpha_strat*R0_func(base,1,i,n,pbt,c,V) + (1-
                                   alpha_strat)*R0_func(base,2,i,n,pbt,c,V)) + rho*base
            pandemic_size.append(1.-fsolve(f,0,args = R)[0])
            
        emnew, fixnew = [], []
        emnew.append(((1-np.sum(mut))*fix(R_wt))) #baseline emergence of wt
        fixnew.append(fix(R_wt)) #baseline fixation of wt
        
        for i in range(n):
            emnew.append(em[i]) #add emergence of other types
            fixnew.append(fixonly[i]) #add fixation of other types
            
        fixnew = np.array(fixnew) #convert to arrays
        emnew = np.array(emnew)
        allem = np.sum(emnew)
        pandemic_size= np.array(pandemic_size)
        sumem = np.zeros((n+1))
        
        for i in range(n+1): #calculation of relative expected pandemic burden
            sumem[n-i] = np.prod(1. - emnew[n-i+1:])**R_IC - np.prod(1. - emnew[n-i:])**R_IC
        
        #final calculations of burden mean, variance and std
        burden = np.sum(sumem*pandemic_size)/pandemic_size[-1]/(1. -(1-fix(base))**R_IC)
        var = np.sum(sumem*pandemic_size**2.)/pandemic_size[-1]**2./(1. -(1-fix(base))**R_IC)**2. - burden**2.
        std = np.sqrt(var)
        
        #record results for this value of rho and alpha
        alpha_result_list.append(burden)
        std_result_list.append(var)
        
    #record relevant results in aggregate lisst
    minalphalist.append(alphaspace[np.array(alpha_result_list).argmin()])
    stdbestlist.append(std_result_list[np.array(alpha_result_list).argmin()])
    mosaiclist.append(alpha_result_list[-1])
    stdmosaiclist.append(std_result_list[-1])
    pyramidlist.append(alpha_result_list[0])
    stdpyramidlist.append(std_result_list[0])
    bestlist.append(min(alpha_result_list))

'''
'''
#creates the data for the inlay in figure (a)
    
alphaspace = np.linspace(0,1,300)
Vspace = np.linspace(0,1,100)
xlabel = r'$\rho$'
minalphalist2 = []
mosaiclist2 = []
pyramidlist2 = []
bestlist2 = []

for rho in rhospace:
    alpha_result_list2 = []
    for alpha_strat in alphaspace:
        alpha_result_list2.append(fix_all_strats(wibes, n, base, ol, 0, c, V, dur, p, q, rho, alpha_strat,N)[0])
        
    minalphalist2.append(alphaspace[np.array(alpha_result_list2).argmin()])
    mosaiclist2.append(alpha_result_list2[-1])
    pyramidlist2.append(alpha_result_list2[0])
    bestlist2.append(min(alpha_result_list2))
 
'''
'''
#plotting results
    
fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize=(15, 5))
newax = fig.add_axes([0.35, 0.3, 0.11, 0.15], anchor='SE')

ax = axs[0]
ax.set_title('(a)',size=20,loc='left')
ax.set_ylabel(r'$p_{est}$', size=20)
ax.set_xlabel(xlabel,size=20)
plt.ylim(-0.05,1.05)
ax.plot(rhospace,bestlist2,'k-', label='optimal: ' + r'$\alpha = \alpha^*$', linewidth = 5)
ax.plot(rhospace,mosaiclist2,'r--', label='1-epitope: ' + r'$\alpha = 1$', linewidth = 5)
ax.plot(rhospace,pyramidlist2,'b--', label='2-epitope: ' + r'$\alpha = 0$', linewidth = 5)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_ylim(-0.05,1.05)

ax = newax
ax.set_xlabel(xlabel,size=20)
ax.set_ylabel(r'$\alpha^*$', size=20)
ax.plot(rhospace[1:],minalphalist2[1:],'k-',linewidth = 2)
ax.tick_params(axis='both', which='major', labelsize=14)

ax = axs[1]
ax.set_title('(b)',size=20,loc='left')
ax.set_ylabel(r'$E \ [s/s_0]$', size=20)
ax.set_xlabel(xlabel,size=20)
ax.set_ylim(-0.05,1.05)
ax.plot(rhospace,bestlist,'k-', label='optimal: ' + r'$\alpha = \alpha^*$', linewidth = 5)
ax.plot(rhospace,mosaiclist,'r--', label='1-epitope: ' + r'$\alpha = 1$', linewidth = 5)
ax.plot(rhospace,pyramidlist,'b--', label='2-epitope: ' + r'$\alpha = 0$', linewidth = 5)

legend = ax.legend(prop={'size':12},title='Strategy:')
ax.tick_params(axis='both', which='major', labelsize=14)
plt.setp(legend.get_title(),fontsize=12)
plt.ylim(-0.05,1.05)
plt.savefig('Figure9.jpg')
plt.show()
