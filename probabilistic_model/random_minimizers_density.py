import numpy as np
import time
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import MaxNLocator

################################################ PARAMETERS
m = 8

M = 1000
w = 10

N_simu = 10 ** 3
max_Z = 100 # expected = 2/(w+1)*M = 181 with M=1000 and w=10

################################################

U = 4**m
mu_theo = (w+1)/2

means = np.zeros(max_Z)

count=0
density = []


print('Generating data...')
t=time.time()
while count < N_simu:

    S = np.random.randint(1, U + 1, M)

    Z = []
    flag = False

    j=0
    prev_position = -1
    prev_value = float('Inf')

    for i in range(M-w):

        if prev_position <i :
            #rescan
            new_position = i + np.argmin(S[i:i + w])
            prev_value = min(S[i:i+w])
            Z.append(new_position-prev_position)
            prev_position=new_position
            j+=1
            flag=True
        else:
            # no rescan
            if S[i+w-1]<prev_value:
                prev_value = S[i+w-1]
                Z.append(i+w-1-prev_position)
                prev_position = i + w-1
                j+=1

    density.append(j/M)

    if j>max_Z:
        sys.stdout.write("\rProgression: %i/%i" % (count + 1, N_simu))
        means+= Z[:max_Z]
        count+=1

print('\n... done in %f seconds' % (time.time()-t))

print('Theoretical density:',1/mu_theo)
print('Empirical density:',np.mean(density))
print('Error percentage : ', 100*np.abs(np.mean(density)-1/mu_theo)*mu_theo)

print('-------')
print('Theoretical Z_2 mean:',(3*np.log(2)-1)/2*w + np.log(2)/2 +1/8 - (w-1)/(32*w*w))
print('Empirical Z_2 mean:',means[1]/N_simu)

z = []
obs_mean = []

for i in range(max_Z):
    z.append(i+1)
    obs_mean.append(means[i]/N_simu)

#####################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
fontsize = 30

#####################################

fig,ax = plt.subplots(figsize=(12, 8))

ax.plot(z,obs_mean)
ax.axhline(y=mu_theo,c='C3')
ax.tick_params(axis='both', which='major', labelsize=fontsize)

plt.show()
fig.tight_layout()
fig.savefig('mean.pdf')

fig,ax = plt.subplots(figsize=(12, 8))

errors = [i-mu_theo for i in obs_mean]

print('Mean of errors :',np.mean(errors))

ax.hist(errors, density=True,bins=100)
ax.axvline(x=np.mean(errors),c='C3')
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.tick_params(axis='both', which='major', labelsize=fontsize)

plt.show()
fig.tight_layout()
fig.savefig('deviation_to_mean.pdf')
