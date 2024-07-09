import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np

#Creates a schematic figure of a pathogen being recognized by antibodies. 

def drawsurface(ax,nstick,nab,nmut,vaccset,mutset, annotate):
    colors2  = ['#377eb8','#e41a1c',
    '#4daf4a',
    '#984ea3',
    '#ff7f00',
    '#ffff33',
    '#a65628',
    '#f781bf']
    
    #nab = 2
    #nmut = 1
    #nstick = 4
    #vaccset = np.random.choice(range(nstick),size=nab,replace=False)
    #mutset = np.random.choice(range(nstick),size=nstick-nmut,replace=False)
    agradius = 1.
    maxangle = 90
    upgrade = 4
    body = upgrade*agradius*nstick*360/(np.pi*maxangle)
    stick = 0.2*body
    d = stick*0.2
    attach = 2*stick
    gamma = np.pi*0.3
    circle = plt.Circle((0, 0), radius=body, ec = 'k',fc='w',lw=2.5)
    ax.add_patch(circle)
    
    i = 0
    #plt.annotate('Pathogen Type \n [1,1,1,1]',(0,0),xytext=(-3,0))
    
    for a in np.linspace(maxangle/(nstick+1),maxangle - maxangle/(nstick+1),nstick):
        alpha = a*np.pi/180
        print alpha
        
        line = plt.Line2D((np.cos(alpha)*body, np.cos(alpha)*(body+stick)), 
                              (np.sin(alpha)*body, np.sin(alpha)*(body+stick)), 
                              lw=2.5, linestyle = '--',color= 'black')
        ax.add_line(line)

        arc = Arc((np.cos(alpha)*(body+stick+agradius),np.sin(alpha)*(body+stick+agradius)),
                   agradius*2, agradius*2, angle = 180+alpha*180/np.pi, 
                   theta1 = -130, theta2 = 130,lw=2.5, linestyle = '--')
        ax.add_patch(arc)
        
        if i in mutset:
            line = plt.Line2D((np.cos(alpha)*body, np.cos(alpha)*(body+stick)), 
                              (np.sin(alpha)*body, np.sin(alpha)*(body+stick)), 
                              lw=2.5,color= 'black')
            ax.add_line(line)
            
            ag = plt.Circle((np.cos(alpha)*(body+stick+agradius),
                             np.sin(alpha)*(body+stick+agradius)),
                            radius = agradius,alpha = 0.5,color = colors2[i])
            ax.add_patch(ag)
            
            arc = Arc((np.cos(alpha)*(body+stick+agradius),np.sin(alpha)*(body+stick+agradius)),
                       agradius*2, agradius*2, angle = 180+alpha*180/np.pi, 
                       theta1 = -130, theta2 = 130,lw=2.5)
            ax.add_patch(arc)
        
        x1 = np.cos(alpha)*(body+stick+attach) + np.sin(alpha)*d/2.
        x2 = np.cos(alpha)*(body+stick+attach) - np.sin(alpha)*d/2.
        y1 = np.sin(alpha)*(body+stick+attach) - np.cos(alpha)*d/2.
        y2 = np.sin(alpha)*(body+stick+attach) + np.cos(alpha)*d/2.
     
        if i in vaccset:
            ag = plt.Circle((np.cos(alpha)*(body+stick+attach),
                             np.sin(alpha)*(body+stick+attach)),
                            radius = agradius,alpha = 0.5,color = colors2[i])
            ax.add_patch(ag)
            
            line = plt.Line2D((np.cos(alpha)*(body+stick+attach) + np.sin(alpha)*d/2.,
                               np.cos(alpha)*(body+2*stick+attach) + np.sin(alpha)*d/2.), 
                              (np.sin(alpha)*(body+stick+attach) - np.cos(alpha)*d/2.,
                               np.sin(alpha)*(body+2*stick+attach) - np.cos(alpha)*d/2.), 
                              lw=2.5,color= 'black')
            ax.add_line(line)
            
            line = plt.Line2D((np.cos(alpha)*(body+stick+attach) - np.sin(alpha)*d/2.,
                               np.cos(alpha)*(body+2*stick+attach) - np.sin(alpha)*d/2.), 
                              (np.sin(alpha)*(body+stick+attach) + np.cos(alpha)*d/2.,
                               np.sin(alpha)*(body+2*stick+attach) + np.cos(alpha)*d/2.), 
                              lw=2.5,color= 'black')
            ax.add_line(line)
            
            line = plt.Line2D((x2,
                               x2 -stick*0.5*np.cos(np.pi*0.5 - alpha - gamma)), 
                              (y2,
                               y2 + stick*0.5*np.sin(np.pi*0.5 - alpha - gamma)), 
                              lw=2.5,color= 'black')
            ax.add_line(line)
            
            line = plt.Line2D((x1,
                               x1 +stick*0.5*np.cos(np.pi*0.5 - alpha + gamma)), 
                              (y1,
                               y1 - stick*0.5*np.sin(np.pi*0.5 - alpha + gamma)), 
                              lw=2.5,color= 'black')
            ax.add_line(line)
        i+=1
    if annotate == 1:
        ax.annotate(r'Epitope', xy=(np.cos(alpha)*body, np.sin(alpha)*body), xytext=(body/4, body/3),
            arrowprops=dict(facecolor='black', shrink=0.08),size=20)
        ax.annotate('Immune \n Response', xy=(x1 + agradius,
                                                    y1 - agradius),
                                                    xytext=(body+1, body+1),
            arrowprops=dict(facecolor='black', shrink=0.08),size=20)
    if annotate == 2:
        ax.annotate('Mutation \nChanged\n Epitope', xy=(np.cos(alpha)*body, np.sin(alpha)*body), xytext=(body/4, body/8),
            arrowprops=dict(facecolor='black', shrink=0.08),size=15)
        
    ax.axis('off')
    ax.axis('scaled')
    ax.axis(xmin = 0, xmax = body+2.*stick+attach, 
                 ymin = 0, ymax = body+2.*stick+attach)
    
#%%
#Create Figure Schematic using n different epitopes, m which are vaccinated, and a pathogen with i mutations
n = 2
m = 2
i = 0
fig, ax = plt.subplots(1,m+1, figsize=(15, 15),dpi=300)

for i in range(0,m+1):
    vaccs = [0,1] #np.random.choice(range(n),size=m,replace=False)
    muts = [0,1][:m-i] #np.random.choice(range(n),size=n-i,replace=False)
    drawsurface(ax[i],n,m,i,vaccs,muts,i+1)
    ax[i].set_title(r'$i = $' +str(i),size=30)

plt.savefig('Figure2.jpg')
#%%


