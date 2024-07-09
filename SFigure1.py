import numpy as np
import matplotlib.pyplot as plt
from pyrasaic_tools import *
import matplotlib.colors as colors

    
#Comparing different within body evolution Regimes: 
#Bernoulli Model, Wright Fisher Model, Infinite N Model

#set parameters
dur = 150
dur_normal  = 10
beta = 0.2
base = beta*dur_normal
q = 10**-8
wibes = 2
p = 0.03
dur_normal  = 10
beta = 0.2
base = beta*dur_normal
ol = 0.
V = 0.8
n= 2
m = 2

T = dur

#define empty data arrays which will later be filled
pathogen_composition1, pathogen_composition2, pathogen_composition3, pathogen_composition4 = [], [],[],[]
pathogen_composition1b, pathogen_composition2b, pathogen_composition3b, pathogen_composition4b = [], [],[],[]
pathogen_composition1c, pathogen_composition2c, pathogen_composition3c, pathogen_composition4c = [], [],[],[]

#simulate evolutionary trajectories using different methods
for t in range(0,T):
    pathogen_composition1.append(convoluted_poissons(m,n,p,q,t)[:(n+1)]) #bernoulli model
pathogen_composition2 = calculate_wibe(n,m,q,wibes,ol,t)[1] #infinite model
pathogen_composition3 = calculate_wf_wibe(n,m,q/p,wibes,ol,T,int(p/q),30)[1] #wright fisher model

for t in range(1,T):
    pathogen_composition1b = pathogen_composition1
    pathogen_composition2b.append(np.mean(pathogen_composition2[:t],axis=0))
    pathogen_composition3b.append(np.mean(pathogen_composition3[:t],axis=0))


#########
#Plotting

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 6),
                         sharey='row', sharex='row')
for i in range(n+1):
    if i == 0:
        axes[0].plot(np.array(pathogen_composition1b)[:,i],label = r"Wildtype $i = 0$")
    else:
        axes[0].plot(np.array(pathogen_composition1b)[:,i],label = r'Mutant $i = $' + str(i))
        
axes[0].legend(prop={'size':12})
axes[0].set_title('Bernoulli Model',size=15)
axes[0].set_title('(a)',size=15,loc ='left')
axes[0].set_xlabel(r'$\tau$',size=20)
axes[0].set_ylabel(r'Pathogen Composition $P_n(i,\tau)$',size=15)
axes[0].tick_params(axis='both', which='major', labelsize=15)
axes[1].plot(pathogen_composition2b)
axes[1].set_title('Infinite Model',size=15)
axes[1].set_title('(b)',size=15,loc ='left')
axes[1].set_xlabel(r'$\tau$',size=20)
axes[1].tick_params(axis='both', which='major', labelsize=15)
axes[2].plot(pathogen_composition3b)
axes[2].set_title('Avg. W. F. Model',size=15)
axes[2].set_title('(c)',size=15,loc ='left')
axes[2].set_xlabel(r'$\tau$',size=20)
axes[2].tick_params(axis='both', which='major', labelsize=15)
plt.show()