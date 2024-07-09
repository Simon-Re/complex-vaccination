import numpy as np
import matplotlib.pyplot as plt
from pyrasaic_tools import *

#size of the monte carlo simulation sample
trys = 1000

#create figure (a) R0-tau
#set parameters
halfways = 0.05
dur = 10
dur_normal  = 10
beta = 0.2
base = beta*dur_normal
base = 2.
q = 10**-10
p = 5*10**-2
N = 1000000
ol= 0.
pbt = 0.
c = 0.
V = 1.
n= 4
wibes = 1
grains =50
spacing = grains/10
mode = 'Conv'

#possible parameter ranges
cs = np.linspace(0,1.,grains)
Ns = np.logspace(1,6,grains)
bases = np.linspace(1.0001,5,grains)
ss = np.linspace(0,1.,grains)
ols = np.linspace(0,1,grains)
pbts = np.linspace(0,1,grains)
ps = np.linspace(0,0.2,grains)
qs = np.linspace(-20,-3,grains)
rhos = np.linspace(0,1./base,grains)
ns = np.arange(1,10+1)
durs = np.linspace(1,29,grains)
Vs = np.linspace(0,1.,grains)
wibess = np.linspace(0.05,20.,grains)
ratios = np.linspace(0.,1.,grains)
Vmesh = [V]

val1 = rhos
val2 = Ns
title1 = r'$\rho$'
title2 = r'$N$'

canvas_worst, canvas_best = np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_worst_m_max, canvas_best_m_max=  np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_worst_m_min, canvas_best_m_min=  np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_best_m, canvas_worst_m=  np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_best_m_len = np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas = np.zeros((n,len(val1),len(val2)))
canvas_best_mix, canvas_best_opt = np.zeros((len(val1),len(val2))), np.zeros((len(val1),len(val2)))
        
for i in range(len(val1)):
    print i
    pbt = val1[i]
    Vmesh = [V]
    for j in range(len(val2)):
        N = val2[j]
        problist = np.array(integrate_through_campaign(mode,wibes,n,
                                                       base,ol,pbt,c,Vmesh,dur,
                                                       p,q,N))
        mixed_strat_probs = []
        concentrations = np.ones(n)
        
        #interatre through multiple instances of alpha values drawn from dirichlet
        for trial in range(trys):
            alphas = np.random.dirichlet(concentrations) #alpha normalized
            mixed_strat_probs.append(fix_all_strats_hd(wibes, n, base, ol, 0, c, V, dur, p, q, pbt, alphas,N)[0])

        canvas_best_m[i,j] = np.min(np.where(problist == problist.min())) +1
        canvas_worst_m[i,j] = np.min(np.where(problist == problist.min())) +1
        canvas_best[i,j] = min(problist)
        canvas_best_mix[i,j] = min(mixed_strat_probs)
        canvas_worst[i,j] = max(problist)
        
        for m in range(n):
            canvas[m,i,j] = problist[m]       
            
        halfways_rel = min(problist) + halfways
        problist = np.heaviside(problist-halfways_rel,1)*problist
        canvas_best_m_min[i,j] = np.min(np.where(problist == problist.min())) +1
        canvas_best_m_max[i,j] = np.max(np.where(problist == problist.min())) +1
        canvas_worst_m_min[i,j] = np.min(np.where(problist == problist.max())) +1
        canvas_worst_m_max[i,j] = np.max(np.where(problist == problist.max())) +1

        if canvas_best[i,j] > halfways:
            canvas_best_m_min[i,j], canvas_worst_m_min[i,j] = np.nan, np.nan
            canvas_best_m_max[i,j], canvas_worst_m_max[i,j] = np.nan, np.nan

#########
#Plotting

fig, axs = plt.subplots(1,2,figsize=(15, 5))


logs = np.log10(np.logspace(np.log10(min(val2)),np.log10(max(val2)),int(np.log10(max(val2)))))
stringlist = []

for i in range(int(max(logs))):
    stringlist.append(r'$10^'+str(int(logs[i])) + '$')
    
    
axi = axs[1]
pos2 = axi.imshow(canvas_best - canvas_best_mix,cmap = 'Greens',origin = 'lower',interpolation='none',
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto',vmin=0)

cbar = plt.colorbar(pos2,ax = axi,shrink=0.8)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel(r'$p_{est}(m^{*})-p_{est}(\alpha^{*})$', rotation=270, size= 15)
cbar.ax.tick_params(labelsize=14)

axi.set_xticks(np.linspace(min(val2),max(val2),int(max(logs))))
axi.set_xticklabels(stringlist)
axi.set_xlabel(title2,size=20)
axi.set_ylabel(title1,size=20)
axi.tick_params(axis='both', which='major', labelsize=15)

cb = canvas_best
cbm = canvas_best_mix


#run another experiment, Figure (b) - rho, N
#code mostly repeats from above

#parmaters
V = 0.8
halfways = 0.05
dur =10
dur_normal  = 10
beta = 0.2
base = beta*dur_normal
base = 2.
q = 10**-10
p = 2*10**-2
N = 10000
ol= 0.
pbt = 0.001
c = 0.
mode = 'Conv'

n= 4
wibes = 1
grains =40
spacing = grains/10

cs = np.linspace(0,1.,grains)
Ns = np.logspace(1,5,grains)
bases = np.linspace(1.0001,5,grains)
ss = np.linspace(0,1.,grains)
ols = np.linspace(0,1,grains)
pbts = np.linspace(0,1,grains)
ps = np.linspace(0,0.2,grains)
qs = np.linspace(-20,-3,grains)
rhos = np.linspace(0,1./base,grains)
ns = np.arange(1,10+1)
durs = np.linspace(10,300,grains)
Vs = np.linspace(0,1.,grains)
wibess = np.linspace(0.05,20.,grains)
ratios = np.linspace(0.,1.,grains)
Vmesh = [V]

val1 = bases
val2 = durs
title1 = r'$R_0$'
title2 = r'$\tau$'

canvas_worst, canvas_best = np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_worst_m_max, canvas_best_m_max=  np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_worst_m_min, canvas_best_m_min=  np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_best_m, canvas_worst_m=  np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_best_m_len = np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas = np.zeros((n,len(val1),len(val2)))
canvas_best_mix, canvas_best_opt = np.zeros((len(val1),len(val2))), np.zeros((len(val1),len(val2)))
        
for i in range(len(val1)):
    print i
    base = val1[i]
    Vmesh = [V]
    for j in range(len(val2)):
        dur = val2[j]
        problist = np.array(integrate_through_campaign(mode,wibes,n,
                                                       base,ol,pbt,c,Vmesh,dur,
                                                       p,q,N))
        mixed_strat_probs = []
        concentrations = np.ones(n)/1.
        for trial in range(trys):
            alphas = np.random.dirichlet(concentrations)
            mixed_strat_probs.append(fix_all_strats_hd(wibes, n, base, ol, 0, c, V, dur, p, q, pbt, alphas,N)[0])

        canvas_best_m[i,j] = np.min(np.where(problist == problist.min())) +1
        canvas_worst_m[i,j] = np.min(np.where(problist == problist.min())) +1
        canvas_best[i,j] = min(problist)
        canvas_best_mix[i,j] = min(mixed_strat_probs)
        canvas_worst[i,j] = max(problist)
        
        for m in range(n):
            canvas[m,i,j] = problist[m]       
            
        halfways_rel = min(problist) + halfways
        problist = np.heaviside(problist-halfways_rel,1)*problist
        canvas_best_m_min[i,j] = np.min(np.where(problist == problist.min())) +1
        canvas_best_m_max[i,j] = np.max(np.where(problist == problist.min())) +1
        canvas_worst_m_min[i,j] = np.min(np.where(problist == problist.max())) +1
        canvas_worst_m_max[i,j] = np.max(np.where(problist == problist.max())) +1
        if canvas_best[i,j] > halfways:
            canvas_best_m_min[i,j], canvas_worst_m_min[i,j] = np.nan, np.nan
            canvas_best_m_max[i,j], canvas_worst_m_max[i,j] = np.nan, np.nan

#########
#Plotting
            
axi = axs[0]
pos2 = axi.imshow(canvas_best - canvas_best_mix,cmap = 'Greens',origin = 'lower',interpolation='none',
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto',vmin=0)
cbar = plt.colorbar(pos2,ax = axi,shrink=0.8)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel(r'$p_{est}(m^{*})-p_{est}(\alpha^{*})$', rotation=270, size= 15)
cbar.ax.tick_params(labelsize=14)

axi.set_xticks(np.linspace(min(val2),max(val2),int(max(logs))))
axi.set_xticklabels(stringlist)
axi.set_xlabel(title2,size=20)
axi.set_ylabel(title1,size=20)
axi.tick_params(axis='both', which='major', labelsize=15)

axs[0].set_title('(a)',size=20,loc='left')
axs[1].set_title('(b)',size=20,loc='left')

plt.show()