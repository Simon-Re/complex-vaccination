from math import factorial
import numpy as np
import matplotlib.pyplot as plt
import itertools
 
#This code creates spherical shapes, which characterize different possible immune types, 
#depending on vaccine administration with n epitopes

def findsubsets(s, n):
    return list(itertools.combinations(s, n))

def nCr(n,r):
    f = factorial
    if n >= r and r >=0:
        return f(n) // f(r) // f(n-r)
    else:
        return 0

#different colorchoices
colors = [
    '#33691E', '#9768D1', '#BDC6FF', '#FC7D49', '#CD0402', '#BEDB39',
    '#FFF176', '#FA5B0F', '#FE4D98', '#4A148C', '#8C8C8C', '#A5FFD2'
]

colors2  = ['#e41a1c',
'#377eb8',
'#4daf4a',
'#984ea3',
'#ff7f00',
'#ffff33',
'#a65628',
'#f781bf']
plt.style.use('default')

n = 4
m = 1

for m in range(1,n+1):
    
    N = nCr(n,m) 
    rows = int(np.sqrt(N))
    piece = 100./m
    f = plt.figure(dpi=300)
    f.set_figwidth(np.ceil(float(N)/rows))
    f.set_figheight(rows)
    pyramids = findsubsets(range(1,n+1),m)
    
    for i in range(1,N+1):
        choicelist= np.zeros(n)
        for j in range(1,n+1):
            if j in pyramids[i-1]:
                choicelist[j-1] = 1.
        y = np.array(piece*choicelist)
        explovec = np.ones(n)*0.1
        ax = plt.subplot( rows, np.ceil(float(N)/rows),i)
        ax.pie(y,colors = colors2, explode = explovec,startangle = 45)
    f.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.savefig('vacctypes' + str(m) +'.jpg',transparent=True)
    plt.show() 



