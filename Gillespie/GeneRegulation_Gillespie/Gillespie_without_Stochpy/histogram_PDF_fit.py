"""
Make a histogram of normally distributed random numbers and plot the
analytic PDF over it
"""
import numpy as np
from scipy.misc import factorial
import matplotlib.pyplot as plt
# import matplotlib.mlab as mlab

# cutoff the first reacmax/max_trajectory reactions for the histogram
all_traj_equi_n = np.zeros((reacmax-reacmax/max_trajectory)*max_trajectory)
all_traj_equi_m = np.zeros((reacmax-reacmax/max_trajectory)*max_trajectory)
j=0
for element in tract_storage:
	k=j+(max_trajectory-1)
	all_traj_equi_n[(j*reacmax/max_trajectory):(k*reacmax/max_trajectory)] = element[1][reacmax/max_trajectory:]
	all_traj_equi_m[(j*reacmax/max_trajectory):(k*reacmax/max_trajectory)] = element[2][reacmax/max_trajectory:]
	j=k

fig = plt.figure()
ax = fig.add_subplot(111)

# the histogram of the data
# draw histogram using the data from all trajectories
plt.title('histogram found from many reactions at steady state and many trajectories')
plt.xlabel('number of proteins of input species n')
plt.ylabel('histogram for the marginal distribution p(n)')
plt.hist(all_traj_equi_n, bins=nmax, range=(0,nmax))

# hist uses np.histogram under the hood to create 'n' and 'bins'.
# np.histogram returns the bin edges, so there will be nmax probability
# density values in n, 51 bin edges in bins and 50 patches.  To get
# everything lined up, we'll compute the bin centers
bincenters = 0.5*(bins[1:]+bins[:-1])
# add a 'best fit' line for the normal PDF
ideal_poisson = np.exp(-g)*1./
l = ax.plot(bincenters, y, 'r--', linewidth=1)

ax.set_xlabel('Smarts')
ax.set_ylabel('Probability')
#ax.set_title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
ax.set_xlim(40, 160)
ax.set_ylim(0, 0.03)
ax.grid(True)

plt.show()