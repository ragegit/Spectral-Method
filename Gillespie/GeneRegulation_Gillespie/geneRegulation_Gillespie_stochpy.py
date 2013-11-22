import stochpy
import stochpy.modules.Analysis as Analysis
import stochpy.modules.PyscesMiniModel as PMM
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pylab as p

smod = stochpy.SSA() # simulation object
smod.Method('FirstReactionMethod') # method with which simulation is done: First Reaction = Gillespie
smod.Model(File='AutoregAndReg.psc',dir='/home/rage/Projects/Python/SpectralMethod_vs_Gillespie/Gillespie/GeneRegulation_Gillespie/') # specific model with which the simulation runs

smod.DoStochSim(method = 'FirstReactionMethod', mode = 'time', end = 500000, trajectories = 1,IsTrackPropensities=False) # stochastic simulation is run. data of 100 trajectories is produced and each is run until time 1000 is reached. time is scaled by rates
#stochpy.SaveInteractiveSession("smod_end1000_traject1000_GeneReg.py")

smod.PlotAverageSpeciesDistributions('n') # plots the distribution of species n averaging over all trajectories
smod.PlotAverageSpeciesDistributions('m')

#use this to extract original simulation data:
num_traject = 1
alldistributions = [] # list will contain the distribution/normalized histogram and its binning of both particles and of each trajectory
for i in xrange(1,num_traject+1):
	smod.GetTrajectoryData(i)
	distributions_i = smod.data_stochsim.species_distributions # returns distributions and bins of both particle types in a tuple, the distribution and its corresponding binning 
	alldistributions.append(distributions_i)
	
## data without a grid:
#num_traject = 1000
#alldata_n = []
#alldata_m = []
#alldistributions = [] # list will contain the distribution/normalized histogram and its binning of both particles and of each trajectory
#for i in xrange(1,num_traject+1):
	#smod.GetTrajectoryData(i)
	#data_m_i = smod.data_stochsim.getSimData('m') # returns an array with reaction times and particle numbers after that specific reaction.
	#data_n_i = smod.data_stochsim.getSimData('n') # = array([[t0,n0],[t1,n1],[t2,n2],...,[ti,ni]])
	#distributions_i = smod.data_stochsim.species_distributions # returns distributions and bins of both particle types in a tuple, the distribution and its corresponding binning 
	#alldata_n.append(data_n_i)
	#alldata_m.append(data_m_i)
	#alldistributions.append(distributions_i)

# prepare for averaging over all trajectories
bin_min_n = np.int(np.min([alldistributions[i][0][0][-1] for i in xrange(num_traject)]))
bin_min_m = np.int(np.min([alldistributions[i][1][0][-1] for i in xrange(num_traject)]))
dist_n = (np.arange(bin_min_n+1), np.mean([alldistributions[i][0][1][:bin_min_n+1] for i in xrange(num_traject)],axis=0)) # average distribution of species n. dist_n = (bins, distribution of n)
dist_m = (np.arange(bin_min_m+1), np.mean([alldistributions[i][1][1][:bin_min_m+1] for i in xrange(num_traject)],axis=0)) # average distribution of species m. dist_m = (bins, distribution of m)
std_n = (np.arange(bin_min_n+1), np.std([alldistributions[i][0][1][:bin_min_n+1] for i in xrange(num_traject)],axis=0)) # standard deviation of the distribution of species n. dist_n = (bins, stds of n)
std_m = (np.arange(bin_min_m+1), np.std([alldistributions[i][1][1][:bin_min_m+1] for i in xrange(num_traject)],axis=0)) # standard deviation of the distribution of species m. dist_m = (bins, stds of m)

plt.figure()
plt.title("distribution of species n from "+str(num_traject)+" trajectories with hillfunction g(n)=q(n)")
#plt.plot(dist_n[1],label="distribution of species n from "+str(num_traject)+" trajectories")
plt.errorbar(dist_n[0],dist_n[1],yerr=1.96/np.sqrt(num_traject)*std_n[1])
#plt.legend()

plt.figure()
#plt.plot(dist_m[1],label="distribution of species m from "+str(num_traject)+" trajectories")
plt.title("distribution of species m from "+str(num_traject)+" trajectories with hillfunction g(n)=q(n)")
plt.errorbar(dist_m[0],dist_m[1],yerr=1.96/np.sqrt(num_traject)*std_m[1])
#plt.legend()
plt.show()

#--------------------------------------find joint probability distribution ... work in progress TODO

smod.GetRegularGrid()

allhist = []#np.zeros((num_traject,bin_min_n,bin_min_m))
for i in xrange(num_traject):
	#smod.GetTrajectoryData(i)
	#smod.GetRegularGrid()
	x = smod.data_stochsim_grid.species[0][i][11:]
	y = smod.data_stochsim_grid.species[1][i][11:]
	bin_min = min(np.max(x),np.max(y))
	hist_n_m_i = np.histogram2d(x, y, bins=bin_min, range=[[0,bin_min],[0,bin_min]], normed=False, weights=None)
	allhist.append(hist_n_m_i) # = np.array([hist_n_m[:bin_min_n][:bin_min_m][1],hist_n_m[:bin_min_n][:bin_min_m][2]])

new_allhisto = [allhist[i][0] for i in np.arange(len(allhist))]
p_mat = np.mean(new_allhisto,axis=0)
x = np.arange(np.shape(p_mat)[0])
X, Y = p.meshgrid(x, x)
fig=plt.figure()
ax = Axes3D(fig)
ax.plot_wireframe(X,Y,p_mat)
ax.set_xlabel('n')
ax.set_ylabel('m')
ax.set_zlabel('p(n,m)')
plt.show()

np.save('Data/initial_hill2/p_n_m_1trajectory_endtime500000.npy',p_mat)
#np.save('Data/initial_poisson_g7/p_n_m_1trajectory_endtime500000.npy',p_mat)

allhist_m = []
plt.figure()
for i in xrange(0,num_traject):
	hist_m_i = np.histogram(smod.data_stochsim_grid.species[1][i][11:],bins=np.max(smod.data_stochsim_grid.species[1][i][11:]),range=(0,np.max(smod.data_stochsim_grid.species[1][i][11:])),normed=True)
	allhist_m.append(hist_m_i)
	plt.plot(hist_m_i[0],label="distribution of species m from the "+str(i)+"th trajectory")
plt.legend()

##----------------------------------------

#mean_n = 0.
#plt.figure()
#for i in xrange(0,num_traject):
	#mean_n += np.mean(alldata_n[i].T[1][20000:])
	#hist_i_n = np.histogram(alldata_n[i].T[1][20000:],bins=np.max(alldata_n[i].T[1][20000:]),range=(0,np.max(alldata_n[i].T[1][20000:])),normed=True)
	#plt.plot(hist_i_n[0],label="distribution of species n from the "+str(i)+"th trajectory")
#plt.legend()
#mean_n = mean_n/num_traject

#mean_m = 0.
#plt.figure()
#for i in xrange(0,num_traject):
	#mean_m += np.mean(alldata_m[i].T[1][20000:])
	#hist_i_m = np.histogram(alldata_m[i].T[1][20000:],bins=np.max(alldata_m[i].T[1][20000:]),range=(0,np.max(alldata_m[i].T[1][20000:])),normed=True)
	#plt.plot(hist_i_m[0],label="distribution of species m from the "+str(i)+"th trajectory")
#plt.legend()
#mean_m = mean_m/num_traject

def find_histogram(species_num,trajectory,starttime,endtime,smod):
	'''find histogram of specific trajectory of the stochpy simulation object smod in a certain timeintervall'''
	smod.GetTrajectoryData(trajectory)
	data = smod.data_stochsim.getDataInTimeInterval(starttime+0.5*(endtime-starttime),0.5*(endtime-starttime)) # data at each reaction time in an array in form: [[t0,n0,m0],[t1,n1,m1],...]
	hist = np.histogram(data.T[species_num],bins=np.max(data.T[species_num]),range=(0,np.max(data.T[species_num])),normed=True)
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

def plot_Gill_vs_Spec(trajectories,starttime,endtime,smod,p_n_Spec,p_m_Spec):
	plt.figure()
	plt.title('marginal probability density of input proteins n')
	plt.xlabel('n')
	plt.ylabel('p(n)')
	plt.plot(plot_all_hist(1,trajectories,starttime,endtime,smod,'Gillespie'))
	plt.plot(p_n_Spec,'r',label='Spectral Method')
	plt.legend()
	plt.show()
	plt.figure()
	plt.title('marginal probability density of output protein m')
	plt.xlabel('m')
	plt.ylabel('p(m)')
	plt.plot(plot_all_hist(2,trajectories,starttime,endtime,smod,'Gillespie'))
	plt.plot(p_m_Spec,'r',label='Spectral Method')
	plt.legend()
	plt.show()

#smod.PrintSpeciesMeans()
#smod.PlotAverageSpeciesDistributions('n')
#smod.PlotAverageSpeciesDistributions('m')

#smod.data_stochsim.getSimData('n')

#smod.GetTrajectoryData(n=10)
#smod.data_stochsim.getDataInTimeInterval(1.,0.5)
#smod.GetTrajectoryData(n=50)
#smod.data_stochsim.getDataInTimeInterval(1.,0.5)
#smod.GetTrajectoryData(n=100)
#smod.data_stochsim.getDataInTimeInterval(1.,0.5)

#p_n_Gill = smod.data_stochsim.species_distributions[0]
#p_m_Gill = smod.data_stochsim.species_distributions[1]

