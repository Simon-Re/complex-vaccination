import numpy as np
import matplotlib.pyplot as plt
from pyrasaic_tools import herd_immunity, R0_func, fix, convoluted_poissons,fix_all_strats

#defines a checkerboard of red and blue tiles
def checkerboard(w, h) :
      re = np.r_[ w*[0,1] ]              # even-numbered rows
      ro = np.r_[ w*[1,0] ]              # odd-numbered rows
      return np.row_stack(h*(re, ro))

#set parameters  
n =2
wibes = 0.1
V = 0.8
dur = 100
q = 10**-10
p = 2*10**-2
N = 10000
rho = 0.001
c = 0
base =1.
ol = 0

#initialize
grains1, grains2 = 100,100 #resution
Vspace = np.linspace(0,1,100)
alphaspace = np.linspace(0,1,100)
basespace = np.linspace(1.1,5.,grains1)
durspace = np.linspace(1,200,grains2)
rhospace = np.linspace(0,1,grains2)
res_alphas, res_bests, res_mosaics, res_pyramids = np.zeros((4,grains1,grains2))

#gather data for figure by iterating through all combination of R0 [base] and tau [dur]
i, j = 0,0
for base in basespace:
    j =0
    for dur in durspace:
        alpha_result_list = []
        for alpha_strat in alphaspace: #different values of alpha
            alpha_result_list.append(fix_all_strats(wibes, n, base, ol, rho, c, V, dur, p, q, rho, alpha_strat,N)[0])
        #optimal alpha
        res_alphas[i,j] = alphaspace[np.array(alpha_result_list).argmin()]
        #results of alpha = 1 (mosaic)
        res_mosaics[i,j] = alpha_result_list[-1]
        #results of alpha = 0 (pyramid)
        res_pyramids[i,j] =alpha_result_list[0]
        #minimally achievable pathogen emergence
        res_bests[i,j] = min(alpha_result_list)
        j +=1
    i +=1

#########
#Plotting

val1 = basespace
val2 = durspace

fig, axs = plt.subplots(1,2,figsize=(15,5))

ax = axs[0]
im = ax.imshow(1- res_bests,cmap = 'Greens', origin = 'lower',interpolation='none',
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto', vmin = 0,vmax=1)
cs= ax.contour(res_alphas, origin = 'lower',interpolation='none', linewidths = 2, 
            colors = 'k', linestyles = 'dashed', levels = 4, 
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto')
ax.clabel(cs, inline=True, fmt = '%2.1f', fontsize=14)#,manual=[(20,1),(20,1.8),(20,2.5),(20,3.2)])
ax.set_ylabel(r'$R_0$',size=20)
ax.set_xlabel(r'$\tau$',size=20)
cbar = plt.colorbar(im,ax = ax,ticks = [0,1],shrink=0.9)
cbar.ax.get_yaxis().labelpad = 20
cbar.set_label(r'$1-p_{est}(\alpha^*)$', rotation=270,size=20, labelpad=20)
cbar.ax.tick_params(labelsize=14)
ax.tick_params(axis='both', labelsize=15)
ax.set_title('(a)',size=20,loc = 'left')

ax= axs[1]
one_mat = (res_mosaics - res_bests)*checkerboard(grains1/2,grains2/2)
one_mat[one_mat ==0] = np.nan
translucience = 1
pos1 = ax.imshow(one_mat,cmap = 'Reds',origin = 'lower',interpolation='none',
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto', alpha = translucience,vmin = 0,vmax=1.)
cbar = plt.colorbar(pos1,ax=ax,ticks = [0,1],shrink=0.8)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel(r'$p_{est}(1) - p_{est}(\alpha^{*})$', rotation=270, size= 15)
cbar.ax.tick_params(labelsize=14)

one_mat = (res_pyramids - res_bests)*np.abs(checkerboard(grains1/2,grains2/2)-1)
one_mat[one_mat ==0] = np.nan
pos1 = ax.imshow(one_mat,cmap = 'Blues',origin = 'lower',interpolation='none',
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto', alpha = translucience,vmin = 0,vmax=1.)
cbar = plt.colorbar(pos1,ax=ax,ticks = [0,1],shrink=0.8)
cbar.ax.get_yaxis().labelpad = 20
cbar.ax.set_ylabel(r'$p_{est}(0) - p_{est}(\alpha^*)$', rotation=270, size= 15)
cbar.ax.tick_params(labelsize=14)
ax.tick_params(axis='both', labelsize=15)
ax.set_ylabel(r'$R_0$',size=20)
ax.set_xlabel(r'$\tau$',size=20)
ax.set_title('(b)',size=20,loc = 'left')

cs= ax.contour(res_bests, origin = 'lower',interpolation='none', linewidths = 2, 
            colors = 'k', linestyles = 'dashed', levels = [0.1], 
           extent=[min(val2),max(val2),min(val1),max(val1)], aspect ='auto')
ax.clabel(cs, inline=True, fmt = '%2.1f', fontsize=14)#,manual=[(20,1),(20,1.8),(20,2.5),(20,3.2)])