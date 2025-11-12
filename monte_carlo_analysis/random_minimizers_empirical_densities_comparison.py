import numpy as np
import time
import sys
import matplotlib.pyplot as plt
from matplotlib import rc

################################################ PARAMETERS
m = 8

M = 10**5
w = 10

N_simu = 10**1

################################################

print('k=',w+m-1)

U = 4**m
mu_theo = (w+1)/2

print('Generating data...')
t=time.time()

x = set()
sd = dict()
ddd = dict()

nb_kmer = dict()
nb_min = dict()
nb_positions = dict()


for j in range(N_simu):
    sys.stdout.write("\rProgression: %i/%i" % (j + 1, N_simu))

    positions = []
    kmer_set = set()
    min_set = set()

    S = np.random.randint(1, U + 1, M)

    prev_position = -1

    for i in range(M - w):

        kmer = tuple(S[i:i+w])
        kmer_set.add(kmer)

        new_position = i +np.argmin(kmer)

        if new_position != prev_position:
            prev_position = new_position
            positions.append(new_position)
            min_set.add(min(kmer))

        if (i+1) % 1000 == 0:
            x.add(i+1)

            if i+1 not in sd.keys():
                sd[i+1] = []
                ddd[i+1] = []
                nb_min[i+1] = []
                nb_positions[i+1] = []
                nb_kmer[i+1]= []

            sd[i+1].append(len(positions)/(i+1))
            ddd[i+1].append(len(min_set)/len(kmer_set))

            nb_min[i+1].append(len(min_set))
            nb_positions[i+1].append(len(positions))
            nb_kmer[i+1].append(len(kmer_set))

print('\n... done in %f seconds' % (time.time()-t))

x = sorted(list(x))

y_sd = []
y_ddd = []
yup_sd = []
ylow_sd = []
yup_ddd = []
ylow_ddd = []

y_min = []
y_positions = []
y_kmers = []

for i in x:
    y_sd.append(np.median(sd[i]))

    l,u = np.quantile(sd[i],[0.025,0.975])
    yup_sd.append(u)
    ylow_sd.append(l)

    y_ddd.append(np.median(ddd[i]))

    l,u =np.quantile(ddd[i],[0.025,0.975])
    yup_ddd.append(u)
    ylow_ddd.append(l)
    
    y_min.append(np.median(nb_min[i]))
    y_kmers.append(np.median(nb_kmer[i]))
    y_positions.append(np.median(nb_positions[i]))


#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

#####################################

fig,ax = plt.subplots(figsize=(12, 8))

ax.axhline(y=1/mu_theo,c='C3',label=r'$2/(w+1)$')

ax.plot(x,y_sd,label=r"$d$",color='C0')
ax.fill_between(x,ylow_sd,yup_sd,color='C0',alpha=.1)

ax.plot(x,y_ddd,label=r"$d^\ast$",color='C1')
ax.fill_between(x,ylow_ddd,yup_ddd,color='C1',alpha=.1)

ax.tick_params(axis='both', which='major', labelsize=fontsize)

leg = ax.legend(fontsize=fontsize)
for line in leg.get_lines():
    line.set_linewidth(4.0)


fig.tight_layout()
plt.show()
fig.savefig('deduplicated_vs_standard.pdf')

fig,ax2 = plt.subplots(figsize=(12, 8))

ax2.plot(x,x,label='Number of windows',ls=(0,(5,5)))
ax2.plot(x,y_kmers,label='Number of distinct k-mers',ls=(5,(5,5)))
ax2.plot(x,y_positions,label='Number of distinct positions')
ax2.plot(x,y_min,label='Number of distinct minimizers')

ax2.tick_params(axis='both', which='major', labelsize=fontsize)

leg = ax2.legend(fontsize=fontsize)
for line in leg.get_lines():
    line.set_linewidth(4.0)

fig.tight_layout()
plt.show()
