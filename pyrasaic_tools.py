import numpy as np
import matplotlib.pyplot as plt
from math import factorial
from matplotlib import cm 
from scipy.stats import poisson
from scipy.stats import expon
from scipy.special import gamma, gammainc, gammaincc

#Main Module, load to run Figure generating code

'''
transition_matrix(n,q,w)

task: 
    - computes a matrix of possible transitions between 2**n pathogen variants.
    - transition is allowed only at hamming distance of 1
params: 
    - n: number of epitopes, int
    - q: mutation probability, float
    - w: degree of conservedness, float > 0
        *if high, some pathogen variants mutate less frequently than others
returns:
    - transition matrix, array [2**n,2**n]
'''

def transition_matrix(n,q, w = 0.):
  matrix = [[0] * (2**n) for i in range(2**n)]
  volatility = np.ones(n)
  
  #defines volatility using a power law proportionality with w
  for j in range(n):
      volatility[j] = (j+1)**(-w)
      
  # Iterate through all the pathogen variants and check distance criterium
  for i in range(2**n):
    for j in range(n):
      k = i ^ (1 << j)
      matrix[i][k] = q*volatility[j]/n
      matrix[k][i] = q*volatility[j]/n
    matrix[i][i] = 1.- sum(matrix[i][:])
  return np.array(matrix)

'''
fitnesses_from_mutN(n,m,s,ol,landscape)

task: 
    - computes different possible fitness landscapes resulting from vaccination
    - as a function of total number of mutated epitopes
    - linear is used primarily
params: 
    - n: number of epitopes, int
    - m: valence of administed vaccine, int
    - s: general slope of landscape change per epitope change, float
    - ol: similarity of different epitopes in terms of antibody response
        *if high, different epitopes result in similar AB response
        *float 0 < x <1
    - landscape: general shape of the landscape, string
returns:
    - 1D fitness landscape, array [n]
'''

def fitnesses_from_mutN(n,m,s,ol,landscape = 'linear'):
    f = np.ones(n+1)*(1.+s)
    for i in range(m+1):
        if landscape == 'linear':
            f[i] = 1. + s*i/m
        elif landscape == 'exponential':
            f[i] = 1. + ol**(m-i)*s*np.sign(i)
        elif landscape == 'exponential2':
            f[i] = 1. + ol**(m-i)*s
        elif landscape == 'saturating':
            f[i] = 1. + (1. - (1.-ol)**(i))*s
        elif landscape == 'lowden':
            f[i] = 1. + s*(float(i)/m + (m-i)*(1. - (1.-ol)**i)/m)
        elif landscape == 'lowden2':
            f[i] = 1. + s*(1. - (m-i)*(1.-ol)**i/m)
    return f

'''
fitnesses(n,m,s,ol)

task: 
    - computes the fitness of a pathogen variant with genotype-bitstring B
    - from the number of 1s in the bitstring, where each 1 indicates 
    - evasion of a specific immune response associated with an epitope.
params: 
    - n: number of epitopes, int
    - m: valence of administed vaccine, int
    - s: general slope of landscape change per epitope change, float
    - ol: similarity of different epitopes in terms of antibody response
        *if high, different epitopes result in similar AB response
        *float 0 < x <1
returns:
    - fitness landscape, array [2**n]

'''


def fitnesses(n,m,s,ol):
    f = np.ones(2**n)
    for i in range(2**n):
        factors = sum(np.array(list(bin(i))[2:][-m:]).astype(float))
        f[i] = 1. + s*factors/m #linear ignoring ol
        f[i] = 1. + s*(1. - (m-factors)*(1.-ol)**factors/m) #linear with ol
    return f


'''
calculate_wibe(n,m,q,s,ol,dur)

task: 
    - integrates evolutionary trajectories of within body evolution (wibe)
    - using the infinite model of population genetics
params: 
    - n: number of epitopes, int
    - m: valence of administed vaccine, int
    - q: mutation probability, float
    - s: slope of fitnesslandscape, float
    - ol: similarity of different epitopes, float
    - dur: length of the disease, the simulation, int
returns:
    - final pathogen composition, wibepop, array [n]
    - pathogen composition time tracked, wibepoplist, array [n,dur]
'''


def calculate_wibe(n,m,q,s,ol,dur):
    f = fitnesses(n,m,s,ol) #compute fitness of each variant
    Q = transition_matrix(n,q) #compute transition p between variants
    x = np.zeros(2**n)
    x[0] = 1. #initial condition, all variants are wild type
    wibepoplist = []
    T = dur
    
    for i in range(1,T):
         #main integration
        meanf = np.sum(x*f) 
        x = np.matmul(Q.T,x*f)/meanf #normalize
        
        wibepop = np.zeros(n+1)
        
        for i in range(2**n):
            mutN = sum(np.array(list(bin(i))[2:][:]).astype(int))
            wibepop[mutN] += x[i] 
        wibepoplist.append(wibepop)
    
    return wibepop, wibepoplist

'''
calculate_wf_wibe(n,m,q,s,ol,dur,N,exp_n)

task: 
    - integrates evolutionary trajectories of within body evolution (wibe)
    - using the wright fisher model of population genetics
params: 
    - n: number of epitopes, int
    - m: valence of administed vaccine, int
    - q: mutation probability, float
    - s: slope of fitnesslandscape, float
    - ol: similarity of different epitopes, float
    - dur: length of the disease, the simulation, int
    - N: number of pathogens in simulation per generation, int
    - exp_n: number of simulations to gather statistics, int
returns:
    - final pathogen composition, meanwibepop, array [n]
    - pathogen composition time tracked, meanwibepoplist, array [n,dur]
'''

def calculate_wf_wibe(n,m,q,s,ol,dur,N,exp_n):
    f = fitnesses(n,m,s,ol)
    Q = transition_matrix(n,q)
    T = dur
    meanwibepop = np.zeros(n+1)
    meanwibepoplist = np.zeros((T-1,n+1))
    
    for exp in range(exp_n):
        x = np.zeros(2**n)
        x[0] = N
        xlist = []
        wibepoplist = []
        
        for i in range(1,T):
            meanf = np.sum(x*f)
            p_x = np.matmul(Q.T,x*f)/meanf
            x = np.random.multinomial(N,p_x)
            xlist.append(x)
            
            wibepop = np.zeros(n+1)
            for i in range(2**n):
                mutN = sum(np.array(list(bin(i))[2:][:]).astype(int))
                wibepop[mutN] += x[i]
            wibepoplist.append(wibepop)
        meanwibepop += np.array(wibepoplist[-1])
        meanwibepoplist += np.array(wibepoplist)
        
    return meanwibepop/exp_n/N, meanwibepoplist/exp_n/N

'''
calculate_wf_wibe(n,r)

task: 
    - computes the binomial coefficient, n over r
'''

def nCr(n,r):
    f = factorial
    if n >= r and r >=0:
        return f(n) // f(r) // f(n-r)
    else:
        return 0

'''
naive_topology(s_over,n_vacc)

task: 
    - computes the probability to be infected from a pathogen if vaccinated
        *the less epitopes are covered by a vaccine, the higher the 
         probability to become infected
    - functionality is old, typically s_over = 0
params: 
    - s_over: degree of epitope similarity, float 0 < x < 1
    - n_vacc: number of epitopes, int
returns:
    - z, which modulates the epidemiological transmission process, array [n_vacc]
'''

def naive_topology(s_over,n_vacc):
    z = []
    for nu_order in range(n_vacc):
        z.append(s_over**nu_order)
    return z

'''
cost(r0,i,beta)

task: 
    - computes a cost function for each mutation, which reduces transmissiblity
params: 
    - i: number of mutations, int
    - beta: cost parameter, float
returns:
    - z, which modulates the epidemiological transmission process, array [n_vacc]
'''

def cost(i,beta):
    return (1-beta)**i

'''
R0_func(base, m_multi, i_evas, n_vacc, s_over, c_val, V)

task: 
    - computes the reproductive rate of a pathogen depending on params
params: 
    - base: R0, the reproductive rate of the wild type if V =0, float
    - m_multi: valence of vaccines administered, int
    - i_evas: number of evaded epitopes by mutant, int
    - n_vacc: total number of available epitopes, int
    - s_over ("pbt"): similarity of epitopes, affects infectiousness (old), float 
    - c_val: reduced transmissibility (cost) of mutations (old), float
    - V: fraction of vaccinated population, float 0<x<1
returns:
    - effective reproductive rate, float
'''
    
def R0_func(base, m_multi, i_evas, n_vacc, s_over, c_val, V):
    failure_sum= 0
#    z = naive_topology(s_over,n_vacc)
#    for nu_order in range(0,n_vacc):
#        failure_sum += z[nu_order]*nCr(n_vacc - i_evas,nu_order)*nCr(i_evas,m_multi-nu_order)
    failure_sum = nCr(i_evas,m_multi)
    failures = 1 - failure_sum/float(nCr(n_vacc,m_multi))
    return base*(1.-failures*V)#*cost(i_evas,c_val)


'''
herd_immunity(base, m_multi, i_evas, n_vacc, s_over, c_val)

task: 
    - computes the herd immunity threshold for a pathogen 
params: 
    - base: R0, the reproductive rate of the wild type if V =0, float
    - m_multi: valence of vaccines administered, int
    - i_evas: number of evaded epitopes by mutant, int
    - n_vacc: total number of available epitopes, int
    - s_over: similarity of epitopes, affects infectiousness (old), float 
    - c_val: reduced transmissibility (cost) of mutations, float
returns:
    - herd immunity threshold, float
'''

def herd_immunity(base, m_multi, i_evas, n_vacc, s_over, c_val):
    failure_sum= 0
    z = naive_topology(s_over,n_vacc)
    for nu_order in range(0,n_vacc):
        failure_sum += z[nu_order]*nCr(n_vacc - i_evas,nu_order)*nCr(i_evas,m_multi-nu_order)
    failures = 1 - failure_sum/float(nCr(n_vacc,m_multi))
    if failures == 0.:
        return np.nan
    return (1. - 1./(base*cost(i_evas,c_val)))/failures


'''
mut(alpha,p,dur,i,n,m)

task: 
    - computes evolution of pathogens using a binomial distribution
params: 
    - alpha:
    - p: probability of mutation and selection, float
    - dur: length of disease, int
    - i: number of mutations, int
    - n: number of available epitopes, int
    - m: administed vaccine valence, int n>=m
returns:
    - probability of a variant with i mutations to evolve
'''

def mut(p,dur,i,n,m):
    mut_p =0
    if m > i:
        mut_p = nCr(dur,i)*p**i*(1-p)**(dur-i)
    if m == i: 
        for j in range(m,dur+1):
            mut_p += nCr(dur,j)*p**j*(1-p)**(dur-j)
    return mut_p


'''
mut_alt(p,dur,i,n,m)

Alternative to mut using a poisson process.
Note that a gammafunction is the sum of a multiple poisson process

'''

def mut_alt(p,dur,i,n,m):
    mut_p =0
    
    #variants with i > m can not evolve
    if m > i:
        mut_p = gammainc(i+1, dur*p)/(p*dur)
    if m == i:
        mut_p = 1
        for k in range(m+1):
            mut_p -= gammainc(k+1, dur*p)/(p*dur)
    return mut_p

'''
convoluted_poissons(m,n,p,q,dur)

task: 
    - computes evolution of pathogens using the combined outcome
    - of two poisson processes aimed to model differential 
    - immune responses to different epitopes (old functionality)
params: 
    - m: vaccine valence, int
    - n: number of epitopes, int
    - p: mutation and selection under memory immune response, float
    - q: mutation and selection under adaptive immune response, float
    - dur: legth of disease, int
    
returns:
    - probability of a variant with i mutations to evolve

'''

def convoluted_poissons(m,n,p,q,dur):
    mut_p_1 = np.zeros(n+1)
    mut_p_2 = np.zeros(n+1)
    
    #Selection Driven (Memory Immune Response)
    for i in range(m+1):  
        if m > i:
            mut_p_1[i] = gammainc(i+1,dur*p)/(p*dur)#poisson.pmf(i,p*t*m/float(n))
        if m == i:
            mut_p_1[i] = 1 - np.sum(mut_p_1[:m+1])
    
    #Mutation Driven (Adaptive Immune Response)
    for j in range(n-m+1):
        if n-m > j:
            mut_p_2[j] = gammainc(j+1,dur*q)/(q*dur)#*(n-m)/float(n))
        elif n-m == j:
            mut_p_2[j] = 1 - np.sum(mut_p_2[:m+1])#*(n-m)/float(n))
        
    return np.convolve(mut_p_1,mut_p_2)

'''
fix(R)

task:
Computes the probability of survival for a single pathogen starting
in a single infected individual, with rep. rate R, in a birth death 
process

'''

def fix(R):
    if R >= 1.:
        return 1.-1./R
    else:
        return 0.
    
'''
fix_all(selection, wibes, m, n, base, ol, rho, c, V, dur, p,q,N)

task: 
    - computes the probability of survival as above but for any pathogen
    - arising in the evolutionary process within a single host
params: 
    - selection: type of evolutionary process used, string
        *simple: bernoulli process, only m=i state computed
        *conv: bernoulli process with adaptive immunity
        *drift: mutations arrive and are selected even if i>m
        *infinite: using the infinite model
        *infinitecum: infinite model, but variants arise at all times
        *wfmod: using the Wright Fisher Model
        *wfmodcum: WF model, but variants arise at all times
    - wibes: strength of selection in WF and infinite model, float
    - m: administed vaccine valence, int
    - n: number of epitopes, int
    - base: basic rep rate of the wildtype, float
    - ol: similarity of epitopes (affects mutations), float
    - rho: fraction of diseased population, float
    - c: cost of mutations, float
    - V: vaccine rollout, float
    - dur: legth of disease, int
    - p: mutation and selection under memory immune response, float
    - q: mutation and selection under adaptive immune response, float
    - N: number of pathogens used in WF model
    
returns:
    - a single probability of ANY but at least 1 variant to surivive
    - survival probability for each variant 

'''
    
def fix_all(selection, wibes, m, n, base, ol, rho, c, V, dur, p,q,N):
    total_fixation = []
    #p = 0.5*np.log(1.+wibes)/np.log(1./q)
    R_IC = rho*N
    
    if selection == 'Simple':
        R = (1-rho)*R0_func(base,m,m,n,0,c,V) + rho*base
        total_fixation.append(1 - mut_alt(p,dur,m,n,m)*fix(R))    
            
    elif selection == 'Conv':
        mut_conv = convoluted_poissons(m,n,p,q,dur)
        for i in range(1,n+1):
            R = (1-rho)*R0_func(base,m,i,n,0,c,V) + rho*base
            total_fixation.append(1 - mut_conv[i]*fix(R))

    elif selection == 'Drift':
        for i in range(1,n+1):
            R = (1-rho)*R0_func(base,m,i,n,0,c,V) + rho*base
            total_fixation.append(1 - mut_alt(p,dur,i,n,n)*fix(R))
   
    elif selection == 'Infinite':
        infinite_after = calculate_wibe(n,m,q,wibes,ol,dur)[0]
        for i in range(1,n+1):
            R = (1-rho)*R0_func(base,m,i,n,0,c,V) + rho*base
            total_fixation.append(1 - infinite_after[i]*fix(R))
                    
    elif selection == 'InfiniteCum':
        infinite_after = np.sum(calculate_wibe(n,m,q,wibes,ol,dur)[1][:dur],axis=0)/dur
        for i in range(1,n+1):
            R = (1-rho)*R0_func(base,m,i,n,0,c,V) + rho*base
            total_fixation.append(1 - infinite_after[i]*fix(R))
            
    elif selection == 'WFmod':
        wf_after = calculate_wf_wibe(n,m,q/p,wibes,ol,dur,int(p/q),10)[0]
        for i in range(1,n+1):
            R = (1-rho)*R0_func(base,m,i,n,0,c,V) + rho*base
            total_fixation.append(1 - wf_after[i]*fix(R))
            
    elif selection == 'WFmodCum':
        wf_after = np.sum(calculate_wf_wibe(n,m,q/p,wibes,ol,dur,int(p/q),10)[1][:dur],axis=0)/dur
        for i in range(1,n+1):
            R = (1-rho)*R0_func(base,m,i,n,0,c,V) + rho*base
            total_fixation.append(1 - wf_after[i]*fix(R))
    
    #computing a single probability of any pathogen to survive
    single_event_fixation = np.prod(total_fixation)
    return 1.- np.exp(R_IC*np.log(single_event_fixation)), total_fixation

'''
fix_list(selection, wibes, n, base, ol, rho, c, V, dur, p,q,N)

task:
similar to fix_all, but iterates through m and creates a list for experimentation
'''

def fix_list(selection, wibes, n, base, ol, rho, c, V, dur, p,q,N):
    fix = []
    variants = []
    for m in range(1,n+1):
        tp, ap = fix_all(selection, wibes, m, n, base, ol, rho, c, V, dur, p,q,N)
        fix.append(tp)
        variants.append(np.sum((1.-np.array(ap))*np.arange(1,n+1))/np.sum(1-np.array(ap)))
    return fix, variants

'''
fix_list(selection, wibes, n, base, ol, rho, c, V, dur, p,q,N)

task:
similar to fix_all, but computes the cummulative probability of survival througout
the whole vaccination campaign

params:
    Vmesh: a list of vaccination rollouts, at which surivival is evaluated, array [*]
'''

def integrate_through_campaign(selection,wibes,n,base,ol, rho, c, Vmesh, dur, p,q,N):
    Vp = 1- np.array(fix_list(selection, wibes, n, base, ol, rho, c, Vmesh[0], dur, p,q,N)[0])
    for V in Vmesh[1:]:
        Vp *= 1-np.array(fix_list(selection,wibes,n,base, ol, rho, c,V,dur,p,q,N)[0]/len(Vmesh[1:]))
    return 1- Vp

'''
fix_all_strats(wibes, n, base, ol, pbt, c, V, dur,p, q, rho, alpha_strat,N)

task:
similar to fix_all, but allows for vaccination strategies at n=2 with varying 
proportions of mosaics and 2-epitope pyramids

params:
    alpha_strat: a parameter that interpolates between m=1 and m=2, float 1 > x > 0

returns:
    list of probability of emergence for each variant, array [n]
    probability that any variant emerges, float
    list of probability of a variant evolves in a host for each variant, array [n]
    probability that any variant evolves, float
    list of probability of a variant, given it evolved, survives in population, array[n]
    probability of survival for any variant, float
'''

def fix_all_strats(wibes, n, base, ol, pbt, c, V, dur,p, q, rho, alpha_strat,N):
    total_emergence = []
    just_mut = []
    just_fix = []
    pbt = 0
    R_IC = rho*N
    n =2
    
    mut_conv = np.zeros((n+1,2*n+1))
    
    #precompute probabilities of mutation an within host evolution
    for m in range(1,n+1):
        mut_conv[m] = convoluted_poissons(m,n,p,q,dur)
     
    #THIS ONLY WORKS FOR n = 2:
    #compute probability of survival on population level
    for i in range(1,n+1):
        R = (1- rho)*(alpha_strat*R0_func(base,1,i,n,pbt,c,V) + (1-
                               alpha_strat)*R0_func(base,2,i,n,pbt,c,V)) + rho*base
        total_emergence.append(1 - alpha_strat*V*mut_conv[1][i]*fix(R) - (1-alpha_strat)*V*mut_conv[2][i]*fix(R))
        just_mut.append(1-alpha_strat*V*mut_conv[1][i] - V*(1-alpha_strat)*mut_conv[2][i])
        just_fix.append(1 -fix(R))
        
    #calculating the probability that at least one variant emerges
    single_event_emergence = np.prod(total_emergence)
    single_event_mutation = np.prod(just_mut)
    single_event_fixation = np.prod(just_fix)
    
    results = 1.- np.exp(R_IC*np.log(single_event_emergence)), 1. - np.array(total_emergence), \
              1.- np.exp(R_IC*np.log(single_event_mutation)), 1. - np.array(just_mut), \
              1.- np.exp(R_IC*np.log(single_event_fixation)), 1. - np.array(just_fix), 

    return results

'''
fix_all_strats_hd(wibes, n, base, ol, pbt, c, V, dur,p, q, rho, alpha_strat,N)

task:
similar to fix_all_strats, but applicable to n > 2 too

params:
    alpha_strat: strategy vector, array [n], normalized
'''

def fix_all_strats_hd(wibes, n, base, ol, pbt, c, V, dur,p, q, rho, alpha_strat,N):
    
    total_emergence = []
    just_mut = []
    just_fix = []
    pbt = 0
    R_IC = rho*N
    
    mut_conv = np.zeros((n+1,2*n+1))
    
    #precompute probabilities of mutation an within host evolution
    for m in range(1,n+1):
        mut_conv[m] = convoluted_poissons(m,n,p,q,dur)

    #WORKS FOR ANY n
    #compute probability of survival on population level
    for i in range(1,n+1):
        R = 0
        for m in range(1,n+1):
            R += (1- rho)*alpha_strat[m-1]*R0_func(base,m,i,n,pbt,c,V) 
        R += rho*base
        
        P= 0
        for m in range(1,n+1):
            P += alpha_strat[m-1]*V*mut_conv[m][i]
            
        total_emergence.append(1 - P*fix(R)) #within host evolution +survival
        just_mut.append(1-P) #within host evolution only
        just_fix.append(1 -fix(R)) #survival only
        
    #calculating the probability that at least one variant emerges
    single_event_emergence = np.prod(total_emergence)
    single_event_mutation = np.prod(just_mut)
    single_event_fixation = np.prod(just_fix)
    
    results = 1.- np.exp(R_IC*np.log(single_event_emergence)), 1. - np.array(total_emergence), \
              1.- np.exp(R_IC*np.log(single_event_mutation)), 1. - np.array(just_mut), \
              1.- np.exp(R_IC*np.log(single_event_fixation)), 1. - np.array(just_fix), 

    return results