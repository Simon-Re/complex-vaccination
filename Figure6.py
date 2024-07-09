import numpy as np
import matplotlib.pyplot as plt
from pyrasaic_tools import herd_immunity, R0_func, fix

fig, axs = plt.subplots(2, 2,figsize=(15, 13))

'''
Compute Herd immunity Figure (b)
'''

#define parameters
n = 4
R0 = 2.
ol = 0.
c = 0.0

colorlist = ['r','b','g','y','k']

#compute herd immunity threshold for different vacccine valences
j = 0
for m in range(1,n+1):
    V_herd_list1 = []
    V_herd_list2 = []
    for i in range(0,n):
        if R0 < 1.:
            V_herd = 0. #define as 0 if pathogen goes extinct
        else:
            V_herd = herd_immunity(R0, m, i, n, ol, c)
        if V_herd > 1.:
            V_herd = np.nan #indefinite if immunity is unachievable with vaccs
        if V_herd < 0.:
            V_herd = 0.
        V_herd_list1.append(V_herd)
            
    axs[0,1].plot(range(0,n),V_herd_list1,colorlist[j]+'*',label= str(m), markersize=10.)
    axs[0,1].plot(range(0,n),V_herd_list1,colorlist[j]+'--',linewidth=3.0,alpha = 0.5)
    j+=1

axs[0,1].set_ylim((1-1./R0-0.1),1)
leg = axs[0,1].legend(prop={'size':15})
leg.set_title('m =',prop={'size':15})
axs[0,1].set_xticks(range(0,n))
axs[0,1].set_xlabel(r'Evaded Epitopes $i$',size=20)
axs[0,1].set_ylabel(r'$V_H(m,i)$',size=20)


'''
Compute Effective Reproductive Rate (a)
'''
#parameters
n = 4
R0 = 2.
V = 0.8
ol = 0.
c = 0.0

axs[0,1].axhline(y = 1., linewidth= 2,color = 'black', linestyle='-.',alpha = 0.5)

colorlist = ['r','b','g','y','k']
j = 0
for m in range(n):
    R0_list1 = []
    R0_list2 = []
    for i in range(0,n+1):
        Ri = R0_func(R0, m+1, i, n, ol, c,V)
        R0_list1.append(Ri)
            
    axs[0,0].plot(range(0,n+1),R0_list1,colorlist[j]+'*',label= str(m+1), markersize=10.)
    axs[0,0].plot(range(0,n+1),R0_list1,colorlist[j]+'--',linewidth=3.0,alpha = 0.5)
    j+=1

legend = axs[0,0].legend(prop={'size':15})
legend.set_title('m=',prop={'size':15})
axs[0,0].set_xticks(range(0,n+1))
axs[0,0].set_xlabel(r'Evaded Epitopes $i$',size=20)
axs[0,0].set_ylabel(r'$R_i(m)$',size=20)


'''
Compute Probability of Establishment (c)
'''

#parameters
n = 4
m = int(n/2.)+1

R0 = 4
V = 0.75
ol = 0.5
c = 0.0

colorlist = ['r','b','g','y','k']
j = 0
for V in np.linspace(0.,1.,5):
    p_list1 = []
    for i in range(0,n+1):
        if R0 < 1.:
            p_est = 0.
        else:
            p_est = fix(R0_func(R0, i, i, n, ol, c,V))
        if p_est < 0.:
            p_est = 0.
        p_list1.append(p_est)
            
    axs[1,0].plot(range(0,n+1),p_list1,colorlist[j]+'*',label= str(V), markersize=10.)#,  drawstyle='steps-mid')
    axs[1,0].plot(range(0,n+1),p_list1,colorlist[j]+'--',linewidth=3.0,alpha = 0.5)
    j+=1

axs[1,0].set_ylim(-0.01,(1-1/R0)+0.1)
legend = axs[1,0].legend(prop={'size':15},loc='best', mode = "expand", ncol = 5)
legend.set_title(r'$V =$',prop = {'size':15})
axs[1,0].set_xticks(range(0,n+1))
axs[1,0].set_xlabel('Evaded Epitopes $i$ ($=m$)',size=20)
axs[1,0].set_ylabel(r'$P_{est}(m)$',size=20)

'''
Compute Probability of Establishment as function of n
'''
#parameters
n_max = 20
R0 = 4
V = 0.75
ol = 0.
c = 0.0

colorlist = ['r','b','g','y','k']

j = 0
for mode in ['[(n+1)/2]','1']:
    fixlist = []
    for n in range(1,n_max+1):
        if mode == '[(n+1)/2]':
            m = int((n)/2.) #optimal value of vaccine valence
        if mode == '1':
            m = 1
        fixed = fix(R0_func(R0, m, m, n, ol, c,V))
        fixlist.append(fixed)
            
    axs[1,1].plot(range(1,n_max+1),fixlist,colorlist[j]+'*',markersize=10.,label=mode)
    axs[1,1].plot(range(1,n_max+1),fixlist,colorlist[j]+'--',linewidth=3.0,alpha = 0.5)
    j+=1
    

#plotting
plt.ylim(-0.1,1)
axs[1,1].axhline(y = 1. - 1./(R0*(1.-V)), linewidth= 2,color = 'black', linestyle='--')
axs[1,1].axhline(y = 1. - 1./(R0), linewidth= 2,color = 'black', linestyle='--')
axs[1,1].axhline(y = 1. - 1./(R0*(1-V*(1-ol))), linewidth= 2,color = 'black', linestyle='--')
legend = axs[1,1].legend(prop={'size':15}, ncol = 2)
legend.set_title(r'm = ',prop={'size':15})
axs[1,1].set_xticks(np.arange(0,n_max+1,4).astype(int))
axs[1,1].set_xlabel(r'Number of Epitopes $n$',size=20)
axs[1,1].set_ylabel(r'$P_{est}(m)$',size=20)

axs[1,0].set_title('(c)',size=20,loc = 'left')
axs[0,0].set_title('(a)',size=20,loc = 'left')
axs[0,1].set_title('(b)',size=20,loc = 'left')
axs[1,1].set_title('(d)',size=20,loc = 'left')

axs[1,0].tick_params(axis='both', which='major', labelsize=15)
axs[0,0].tick_params(axis='both', which='major', labelsize=15)
axs[0,1].tick_params(axis='both', which='major', labelsize=15)
axs[1,1].tick_params(axis='both', which='major', labelsize=15)

plt.show()
