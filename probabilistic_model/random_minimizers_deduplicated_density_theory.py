import time
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
from scipy.special import comb

from decimal import *

getcontext().prec = 100

def deduplicated_density(k,m):
    U = 4**m
    V = 4**k
    w = Decimal(k-m+1)

    res = Decimal(0)
    for n in range(1,V+1):
        val= Decimal(0)
        for r in range(1,U+1):
            val += (Decimal(1) - (Decimal(r)/Decimal(U))**w + (Decimal(r-1)/Decimal(U))**w)**Decimal(n)
        res+= (Decimal(U)-val) /Decimal(n) * Decimal(comb(V,n,True))

    return float(res / Decimal(2**V-1))

################################################ PARAMETERS

min_k = 2
max_k = 8

################################################

values = {}
values_bis = {}

print('Beginning computation ...')
t=time.time()

count=0

for k in range(min_k,max_k+1):
    for m in range(1,k+1):
        sys.stdout.write("\rProgression: %i/%i" % (count + 1, max_k*(max_k+1)/2))
        w = k-m+1
        values[(k,w)] = deduplicated_density(k,m)
        count+=1

print('\n... done in %f seconds' % (time.time()-t))

#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

#####################################

ws = list(range(1,max_k+1))
densities=[]
lb = []
for w in ws:
    densities.append(2/(w+1))
    lb.append(1/w)

kdic={}

for k in range(min_k,max_k+1):
    kdic[k]=[]
    for w in range(1,k+1):
        kdic[k].append(values[(k,w)])


fig,ax = plt.subplots(figsize=(12, 8))

for k in range(min_k,max_k+1):
    ax.plot(ws[:len(kdic[k])],kdic[k],marker='o',ls='dotted',label=r'$k='+str(k)+'$')

ax.plot(ws,densities,label=r'$2/(w+1)$')
ax.set_xlabel(r'$w$',fontsize=fontsize)
ax.set_ylabel(r'$d^\ast$',fontsize=fontsize)

ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.tick_params(axis='both', which='major', labelsize=fontsize)
leg=ax.legend(fontsize=fontsize,ncol=2)
for line in leg.get_lines():
    line.set_linewidth(4.0)

plt.show()
fig.tight_layout()
fig.savefig('deduplicated_density_random_minimizers.pdf')
