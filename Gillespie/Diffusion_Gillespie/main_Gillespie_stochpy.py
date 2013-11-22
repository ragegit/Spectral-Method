# one dimensional diffusion with two sites is simulated with Gillespie. The resulting probability distribution after a chosen number of trajectories is compared with the expected Bernoulli distribution.
import stochpy
import numpy as np
import math
import matplotlib.pyplot as plt

smod = stochpy.SSA()
smod.Method('FirstReactionMethod')
smod.Model(File='TwoCellDiffusion.psc',dir='/home/rage/Projects/Python/SpectralMethod_vs_Gillespie/Gillespie/Diffusion_Gillespie/')

smod.DoStochSim(method = 'FirstReactionMethod', mode = 'time', end = 1000,
trajectories = 10,IsTrackPropensities=False)
smod.ShowOverview()
smod.PlotAverageSpeciesDistributions('n1')
smod.PlotAverageSpeciesDistributions('n2')
smod.PlotSpeciesDistributions('n1')
smod.PlotSpeciesDistributions('n2')
smod.PrintSpeciesMeans()
#stochpy.SaveInteractiveSession("interactiveSession1.py")

# plot the analytically calculated distribution against the simulation result
def Bernoulli_dist(k,N):
	"""function to dray the distribution"""
	return 1./2.**N*math.factorial(N)/math.factorial(N-k)*1./math.factorial(k)

def find_histogram(species_num,trajectory,starttime,endtime,smod):
	'''find histogram of specific trajectory of the stochpy simulation object smod in a certain timeintervall. the function returns a tuple of two arrays of the same length. the first array contains all possible particle numbers and the second their average occurences'''
	smod.GetTrajectoryData(trajectory)
	data = smod.data_stochsim.getDataInTimeInterval(starttime+0.5*(endtime-starttime),0.5*(endtime-starttime)) # data at each reaction time in an array in form: [[t0,n0,m0],[t1,n1,m1],...]
	hist = np.histogram(data.T[1],bins=np.max(data.T[species_num]),range=(0,np.max(data.T[species_num])),normed=True)
	return (hist[1][:-1],hist[0])

def hist_list(species_num,trajectories,starttime,endtime,smod):
	''''''
	alltra = []
	for tra in np.arange(1,trajectories):
		hist = find_histogram(species_num,tra,starttime,endtime,smod)[1]
		alltra.append(hist)
	maxnum = np.max(map(np.shape,alltra))
	new_alltra = [np.pad(a, (0,maxnum-np.shape(a)[0]), 'constant') for a in alltra]
	return new_alltra

def plot_all_hist(species_num,trajectories,starttime,endtime,smod,text):
	'''returns distribution of the species specified by its number in the smod.data_stochsim.getDataInTimeInterval(...) array'''
	new_alltra = hist_list(species_num,trajectories,starttime,endtime,smod)
	plt.plot(np.sum(new_alltra,axis=0)/np.sum(new_alltra),label=text)
	return np.sum(new_alltra,axis=0)/np.sum(new_alltra)
	
plt.figure()
plt.title('marginal probability density of input proteins n1')
plt.xlabel('n1')
plt.ylabel('p(n1)')
plt.plot(plot_all_hist(1,trajectories,starttime,endtime,smod,'Gillespie'))
plt.plot(np.arange(101),map(Bernoulli_dist,np.arange(101),np.ones(101)*100.),label='Bernoulli')
plt.legend()
plt.show()

plt.figure()
plt.title('marginal probability density of input proteins n1')
plt.xlabel('n2')
plt.ylabel('p(n2)')
plt.plot(plot_all_hist(2,trajectories,starttime,endtime,smod,'Gillespie'))
plt.plot(np.arange(101),map(Bernoulli_dist,np.arange(101),np.ones(101)*100.),label='Bernoulli')
plt.legend()
plt.show()
