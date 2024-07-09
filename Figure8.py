import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm 
from pyrasaic_tools import *

def checkerboard(w, h,g_size) :
      re = np.r_[ w/g_size*([np.nan]*g_size+ [1]*g_size) ]              # even-numbered rows
      ro = np.r_[ w/g_size*([1]*g_size+ [np.nan]*g_size) ]              # odd-numbered rows
      return np.row_stack(h/g_size*([re]*g_size, [ro]*g_size))
  
def checkerboard2(w, h,g_size) :
      re = np.r_[ w/g_size*([1]*g_size+ [np.nan]*g_size) ]              # even-numbered rows
      ro = np.r_[ w/g_size*([np.nan]*g_size+ [1]*g_size) ] 
      return np.row_stack(h/g_size*([re]*g_size, [ro]*g_size))

#set parameters
halfways = 0.1 #used as an allowed "margin of error" in evaluating optima
dur =10
dur_normal  = 10
beta = 0.2
base = beta*dur_normal
base = 2.
q = 10**-10
p = 2*10**-2
N = 10000
rho = 0.001
ol= 0.
pbt = 0.001
c = 0.
V = 0.8
n= 4
wibes = 1
grains =400
spacing = grains/100
mode = 'Conv'

#initialize different parameter spaces of interest
cs = np.linspace(0,1.,grains)
Ns = np.logspace(1,5,grains)
bases = np.linspace(1.0001,5,grains)
ss = np.linspace(0,1.,grains)
ols = np.linspace(0,1,grains)
pbts = np.linspace(0,1,grains)
ps = np.linspace(0,0.2,grains)
qs = np.linspace(-20,-3,grains)
rhos = np.linspace(0,1.,grains)
ns = np.arange(1,10+1)
durs = np.linspace(10,300,grains)
Vs = np.linspace(0,1.,grains)
wibess = np.linspace(0.05,20.,grains)
ratios = np.linspace(0.,1.,grains)
Vmesh = [V] #can be extended, now when iterating through campaign evalueates only at value V

#define the parameters of interest and their title
val1 = bases
val2 = durs
title1 = r'$R_0$'
title2 = r'$\tau$'

#initialize array which will later be filled with data for plotting
canvas_worst, canvas_best = np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_worst_m_max, canvas_best_m_max=  np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_worst_m_min, canvas_best_m_min=  np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_best_m, canvas_worst_m=  np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas_best_m_len = np.zeros((len(val1),len(val2))),np.zeros((len(val1),len(val2)))
canvas = np.zeros((n,len(val1),len(val2)))
        
for i in range(len(val1)):
    base = val1[i] #change this to switch parameter of interest
    Vmesh = [V]
    for j in range(len(val2)):
        dur = val2[j] #change this to switch parameter of interest
        problist = np.array(integrate_through_campaign(mode,wibes,n,
                                                       base,ol,rho,c,Vmesh,dur,
                                                       p,q,N))
        #collect relevant data from calculation results, such as
        #the vaccine valence that minimized emerngece
        canvas_best_m[i,j] = np.min(np.where(problist == problist.min())) +1
        #the vaccine valence that maximizes emergence
        canvas_worst_m[i,j] = np.min(np.where(problist == problist.min())) +1
        #the lowest possibly achievable emergence
        canvas_best[i,j] = min(problist)
        #the highest possibly achievable emergence
        canvas_worst[i,j] = max(problist)
        
        for m in range(n):
            canvas[m,i,j] = problist[m]       
            
        halfways_rel = min(problist) + halfways #smooth results
        problist = np.heaviside(problist-halfways_rel,1)*problist
        
        #minimal cooptimal valence
        canvas_best_m_min[i,j] = np.min(np.where(problist == problist.min())) +1
        #maximal cooptimal valence
        canvas_best_m_max[i,j] = np.max(np.where(problist == problist.min())) +1
        #minimal anti-cooptimal valence
        canvas_worst_m_min[i,j] = np.min(np.where(problist == problist.max())) +1
        #maximal anti-cooptimal valence
        canvas_worst_m_max[i,j] = np.max(np.where(problist == problist.max())) +1
        
        #if the best emergence is above halfways, campaign is considered as failure
        if canvas_best[i,j] > halfways:
            canvas_best_m_min[i,j], canvas_worst_m_min[i,j] = np.nan, np.nan
            canvas_best_m_max[i,j], canvas_worst_m_max[i,j] = np.nan, np.nan

#######
# Plotting first and second row

fig, axs = plt.subplots(3,2,figsize=(15, 20))

one_mat = (canvas[n-1]-canvas_best)*checkerboard(grains/2,grains/2,spacing)
translucience = 1
pos1 = axs[0,1].imshow(one_mat,cmap = 'Blues',origin = 'lower',interpolation='none',
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto', alpha = translucience,vmin = 0,vmax=1.)
cbar = plt.colorbar(pos1,ax=axs[0,1],ticks = [0,1],shrink=0.8)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel(r'$p_{est}(n) - p_{est}(m^{*})$', rotation=270, size= 15)
cbar.ax.tick_params(labelsize=14)

one_mat = (canvas[0]-canvas_best)*checkerboard2(grains/2,grains/2,spacing)
pos1 = axs[0,1].imshow(one_mat,cmap = 'Reds',origin = 'lower',interpolation='none',
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto', alpha = translucience,vmin = 0,vmax=1.)
cbar = plt.colorbar(pos1,ax=axs[0,1],ticks = [0,1],shrink=0.8)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel(r'$p_{est}(1) - p_{est}(m^{*})$', rotation=270, size= 15)
cbar.ax.tick_params(labelsize=14)

axs[0,1].set_xticks(np.linspace(min(val2),max(val2),2))
axs[0,1].set_xlabel(title2,size=20)
axs[0,1].tick_params(axis='both', which='major', labelsize=15)

pos2 = axs[0,0].imshow(1-canvas_best,cmap = 'Greens',origin = 'lower',interpolation='none',
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto', vmin = 0,vmax=1)
cbar = plt.colorbar(pos2,ax = axs[0,0],ticks=[0,1],shrink=0.8)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel(r'$1-p_{est}(m^{*})$', rotation=270, size= 15)
cbar.ax.tick_params(labelsize=14)

axs[0,0].set_xticks(np.linspace(min(val2),max(val2),2))
axs[0,0].set_xlabel(title2,size=20)
axs[0,0].set_ylabel(title1,size=20)
axs[0,0].tick_params(axis='both', which='major', labelsize=15)

newcolors = cm.get_cmap('RdYlBu_r', n)
fmt = '%2.2f'
cs = axs[1,0].contour(canvas_best,origin = 'lower',interpolation='none',levels= [halfways], 
                 colors='black',linestyles='-',
               linewidths = 2,extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto')
cs = axs[1,0].contour(canvas_best_m_min,origin = 'lower',interpolation='none',levels= n-1, 
                 colors='black',linestyles='-',
               linewidths = 2,extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto')

axs[1,0].text(150,3.5, r'$p_{est}(m^*) > 0.1$', fontsize=15)

pos = axs[1,0].imshow(canvas_best_m_min,cmap = newcolors, alpha = 0.5,origin = 'lower',interpolation='none',
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto', vmin = 1,vmax=n)
cbar = plt.colorbar(pos,ax = axs[1,0],ticks=range(1,n+1),shrink=0.8)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel('Minimal Valence \n Co-Optimal Vaccine', rotation=270, size= 15,labelpad = 40)
cbar.ax.tick_params(labelsize=14)

axs[1,0].set_xticks(np.linspace(min(val2),max(val2),2))
axs[1,0].set_xlabel(title2,size=20)
axs[1,0].set_ylabel(title1,size=20)
axs[1,0].tick_params(axis='both', labelsize=15)

fmt = '%2.2f'
cs = axs[1,1].contour(canvas_best,origin = 'lower',interpolation='none',levels= [halfways], 
                 colors='black',linestyles='-',
               linewidths = 2,extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto')
cs = axs[1,1].contour(canvas_best_m_max,origin = 'lower',interpolation='none',levels= 2, 
                 colors='black',linestyles='-',
               linewidths = 2,extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto')

axs[1,1].text(150,3.5, r'$p_{est}(m^*) > 0.1$', fontsize=15)

pos = axs[1,1].imshow(canvas_best_m_max,cmap = newcolors, alpha = 0.5,origin = 'lower',interpolation='none',
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto', vmin = 1,vmax=n)
cbar = plt.colorbar(pos,ax = axs[1,1],ticks=range(1,n+1),shrink=0.8)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel('Maximal Valence \n Co-Optimal Vaccine' ,rotation=270, size= 15,labelpad = 40)
cbar.ax.tick_params(labelsize=14)

axs[1,1].set_xticks(np.linspace(min(val2),max(val2),2))
axs[1,1].set_xlabel(title2,size=20)
axs[1,1].tick_params(axis='both', labelsize=15)

axs[0,0].set_title('(a)',size=20,loc='left')
axs[0,1].set_title('(b)',size=20,loc='left')
axs[1,0].set_title('(c)',size=20,loc='left')
axs[1,1].set_title('(d)',size=20,loc='left')

#%%

#Run another experiment for figure (e) and (f)

#set parameters
halfways = 0.1
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
V = 0.8
n= 6
wibes = 1
grains =30
mode = 'Conv'

#define parameters spaces once more
cs = np.linspace(0,1.,grains)
Ns = np.logspace(1,5,grains)
bases = np.linspace(0.5,5.,grains)
ss = np.linspace(0,1.,grains)
ols = np.linspace(0,1,grains)
pbts = np.linspace(0,1,grains)
ps = np.linspace(0,0.2,grains)
qs = np.linspace(-20,-3,grains)
rhos = np.linspace(0,1.,grains)
ns = np.arange(1,10+1)
durs = np.linspace(10,300,grains)#.astype(int)
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

for n in range(1,n+1):
    for i in range(len(val1)):
        #print i
        base = val1[i]
        Vmesh = [V]
        for j in range(len(val2)):
            dur = val2[j]
            #problist = np.array(fix_list(mode,n,base,s,V,dur,p))
            problist = np.array(integrate_through_campaign(mode,wibes,n,
                                                           base,ol,pbt,c,Vmesh,dur,
                                                           p,q,N))
            canvas_best_m[i,j] = np.min(np.where(problist == problist.min())) +1
            canvas_worst_m[i,j] = np.min(np.where(problist == problist.min())) +1
            canvas_best[i,j] = min(problist)
            canvas_worst[i,j] = max(problist)
            
            for m in range(n):
                canvas[m,i,j] = problist[m]       
                
            halfways_rel = max(problist)*halfways
            problist = np.heaviside(problist-halfways_rel,1)*problist
            #problist = np.flip(problist)
            canvas_best_m_min[i,j] = np.min(np.where(problist == problist.min())) +1
            canvas_best_m_max[i,j] = np.max(np.where(problist == problist.min())) +1
            canvas_worst_m_min[i,j] = np.min(np.where(problist == problist.max())) +1
            canvas_worst_m_max[i,j] = np.max(np.where(problist == problist.max())) +1
            #canvas_best_m_len[i,j] = np.where(problist == problist.max()).len()
            if canvas_best[i,j] > halfways_rel:
                canvas_best_m_min[i,j], canvas_worst_m_min[i,j] = np.nan, np.nan
                canvas_best_m_max[i,j], canvas_worst_m_max[i,j] = np.nan, np.nan

    #create contour plot (e)
    if n ==4:
        col = 'red'
    elif n == 1:
        col = 'black'
    else:
        col = 'grey'
    cs = axs[2,0].contour(canvas_best,origin = 'lower',interpolation='none',levels= [halfways], 
                 colors=col,linestyles='--',
               linewidths = 2,extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto')
    plt.clabel(cs, inline=True, fmt = 'n = ' + str(n), fontsize=14,manual = [(200,1.5)])
    
    #axs[0,0].ticklabel_format(style='sci',scilimits = (0,0))
    axs[2,0].set_xticks(np.linspace(min(val2),max(val2),2))
    axs[2,0].set_xlabel(title2,size=20)
    axs[2,0].set_ylabel(title1,size=20)
    axs[2,0].tick_params(axis='both', which='major', labelsize=15)
    
n = 4
V = 0.8
for V in [0.,0.4,0.6,0.8,0.9]:
    for i in range(len(val1)):
        #print i
        base = val1[i]
        Vmesh = [V]
        for j in range(len(val2)):
            dur = val2[j]
            #problist = np.array(fix_list(mode,n,base,s,V,dur,p))
            problist = np.array(integrate_through_campaign(mode,wibes,n,
                                                           base,ol,pbt,c,Vmesh,dur,
                                                           p,q,N))
            canvas_best_m[i,j] = np.min(np.where(problist == problist.min())) +1
            canvas_worst_m[i,j] = np.min(np.where(problist == problist.min())) +1
            canvas_best[i,j] = min(problist)
            canvas_worst[i,j] = max(problist)
            
            for m in range(n):
                canvas[m,i,j] = problist[m]       
                
            halfways_rel = max(problist)*halfways
            problist = np.heaviside(problist-halfways_rel,1)*problist
            canvas_best_m_min[i,j] = np.min(np.where(problist == problist.min())) +1
            canvas_best_m_max[i,j] = np.max(np.where(problist == problist.min())) +1
            canvas_worst_m_min[i,j] = np.min(np.where(problist == problist.max())) +1
            canvas_worst_m_max[i,j] = np.max(np.where(problist == problist.max())) +1
            
            if canvas_best[i,j] > halfways_rel:
                canvas_best_m_min[i,j], canvas_worst_m_min[i,j] = np.nan, np.nan
                canvas_best_m_max[i,j], canvas_worst_m_max[i,j] = np.nan, np.nan

    #create contour plot (f) at which point emergence hits the value halfways
    if V == 0.8:
        col = 'red'
    elif V == 0.:
        col = 'black'
    else:
        col = 'gray'
    cs = axs[2,1].contour(canvas_best,origin = 'lower',interpolation='none',levels= [halfways], 
                 colors=col,linestyles='--',
               linewidths = 2,extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto')
 
    if V == 0:
        plt.clabel(cs, inline=True, fmt = r'$V \approx $' + str(V), fontsize=14,manual = [(200.,0)])
    else:
        plt.clabel(cs, inline=True, fmt = r'$V = $' + str(V), fontsize=14,manual = [(200.,0)])
 
    axs[2,1].set_xticks(np.linspace(min(val2),max(val2),2))
    axs[2,1].set_xlabel(title2,size=20)
    axs[2,1].tick_params(axis='both', which='major', labelsize=15)

axs[2,0].set_title('(e)',size=20,loc = 'left')
axs[2,1].set_title('(f)',size=20,loc = 'left')