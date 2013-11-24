# A two species gene regulation cascade takes place at two cells which are coupled via diffusion. The system is simulated with Gillespie.
import stochpy
import numpy as np
import math
import matplotlib.pyplot as plt
import pylab as p
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

smod = stochpy.SSA()
smod.Method('FirstReactionMethod')
smod.Model(File='TwoCellDiffusionGeneRegulation.psc',dir='/home/rage/Projects/Python/SpectralMethod_vs_Gillespie/Gillespie/GeneRegulation_Diffusion_Gillespie/')
smod.DoStochSim(method = 'FirstReactionMethod', mode = 'time', end = 100000,trajectories = 1,IsTrackPropensities=False)
smod.PlotAverageSpeciesDistributions(['n1','n2'])
smod.PlotAverageSpeciesDistributions(['m1','m2'])

#---------------

smod.GetRegularGrid()

allhist = []
for i in xrange(num_traject):
	u = smod.data_stochsim_grid.species[0][i][11:]
	v = smod.data_stochsim_grid.species[1][i][11:]
	w = smod.data_stochsim_grid.species[2][i][11:]
	x = smod.data_stochsim_grid.species[3][i][11:]
	#bin_min_n = min(np.max(x),np.max(u),np.max(v),np.max(w))
	#bin_min_m = min(np.max(x),np.max(u),np.max(v),np.max(w))
	hist_n1_n2_m1_m2_i = np.histogramdd([u,v,w,x], bins=(np.max(u),np.max(v),np.max(w),np.max(x)), normed=True, weights=None)
	allhist.append(hist_n1_n2_m1_m2_i) # = np.array([hist_n_m[:bin_min_n][:bin_min_m][1],hist_n_m[:bin_min_n][:bin_min_m][2]])

np.save('Data/initial_hill2/p_n1_n2_m1_m2_1trajectory_endtime100000.npy',hist_n1_n2_m1_m2_i[0])

#TODO
new_allhisto = [allhist[i][0] for i in np.arange(len(allhist))]
p_mat = np.mean(new_allhisto,axis=0)
x = np.arange(np.shape(p_mat)[0])
X, Y = p.meshgrid(x, x)
fig=plt.figure()
ax = Axes3D(fig)
ax.plot_wireframe(X,Y,p_mat)
ax.set_xlabel('n1')
ax.set_ylabel('n2')
ax.set_zlabel('p(n1,n2)')
plt.show()

np.save('Data/initial_hill2/p_n1_n2_1trajectory_endtime100000.npy',p_mat)
#np.save('Data/initial_poisson_g7/p_n_m_1trajectory_endtime500000.npy',p_mat)




#num_traject = 1000
#alldata_n1 = []
#alldata_m1 = []
#alldata_n2 = []
#alldata_m2 = []
#alldistributions = [] # list will contain the distribution/normalized histogram and its binning of both particles and of each trajectory
#for i in xrange(1,num_traject+1):
	#smod.GetTrajectoryData(i)
	#data_m1_i = smod.data_stochsim.getSimData('m1') # returns an array with reaction times and particle numbers after that specific reaction.
	#data_n1_i = smod.data_stochsim.getSimData('n1') # = array([[t0,n0],[t1,n1],[t2,n2],...,[ti,ni]])
	#data_m2_i = smod.data_stochsim.getSimData('m2') # returns an array with reaction times and particle numbers after that specific reaction.
	#data_n2_i = smod.data_stochsim.getSimData('n2') # = array([[t0,n0],[t1,n1],[t2,n2],...,[ti,ni]])
	#distributions_i = smod.data_stochsim.species_distributions # returns distributions and bins of all particle types in a four-tuple, the distribution and its corresponding binning 
	#alldata_n1.append(data_n1_i)
	#alldata_m1.append(data_m1_i)
	#alldata_n1.append(data_n2_i)
	#alldata_m1.append(data_m2_i)
	#alldistributions.append(distributions_i)

## prepare for averaging over all trajectories
#bin_min_n1 = np.int(np.min([alldistributions[i][0][0][-1] for i in xrange(num_traject)]))
#bin_min_m1 = np.int(np.min([alldistributions[i][1][0][-1] for i in xrange(num_traject)]))
#dist_n1 = (np.arange(bin_min_n1+1), np.mean([alldistributions[i][0][1][:bin_min_n1+1] for i in xrange(num_traject)],axis=0)) # average distribution of species n. dist_n = (bins, distribution of n)
#dist_m1 = (np.arange(bin_min_m1+1), np.mean([alldistributions[i][1][1][:bin_min_m1+1] for i in xrange(num_traject)],axis=0)) # average distribution of species m. dist_m = (bins, distribution of m)
#dist_n2 = (np.arange(bin_min_n2+1), np.mean([alldistributions[i][2][1][:bin_min_n2+1] for i in xrange(num_traject)],axis=0)) # average distribution of species n. dist_n = (bins, distribution of n)
#dist_m2 = (np.arange(bin_min_m2+1), np.mean([alldistributions[i][3][1][:bin_min_m2+1] for i in xrange(num_traject)],axis=0)) # average distribution of species m. dist_m = (bins, distribution of m)
#std_n1 = (np.arange(bin_min_n1+1), np.std([alldistributions[i][0][1][:bin_min_n1+1] for i in xrange(num_traject)],axis=0)) # standard deviation of the distribution of species n1. dist_n1 = (bins, stds of n1)
#std_m1 = (np.arange(bin_min_m1+1), np.std([alldistributions[i][1][1][:bin_min_m1+1] for i in xrange(num_traject)],axis=0)) # standard deviation of the distribution of species m1. dist_m1 = (bins, stds of m1)
#std_n2 = (np.arange(bin_min_n2+1), np.std([alldistributions[i][2][1][:bin_min_n2+1] for i in xrange(num_traject)],axis=0)) # standard deviation of the distribution of species n2. dist_n2 = (bins, stds of n2)
#std_m2 = (np.arange(bin_min_m2+1), np.std([alldistributions[i][3][1][:bin_min_m2+1] for i in xrange(num_traject)],axis=0)) # standard deviation of the distribution of species m2. dist_m2 = (bins, stds of m2)


#plt.figure()
#plt.title("distribution of species n1 from "+str(num_traject)+" trajectories with hillfunction g(n1)=q(n1)")
##plt.plot(dist_n[1],label="distribution of species n from "+str(num_traject)+" trajectories")
#plt.errorbar(dist_n1[0],dist_n1[1],yerr=std_n1[1])
##plt.legend()
#plt.figure()
##plt.plot(dist_m[1],label="distribution of species m from "+str(num_traject)+" trajectories")
#plt.title("distribution of species m from "+str(num_traject)+" trajectories with hillfunction g(n1)=q(n1)")
#plt.errorbar(dist_m1[0],dist_m1[1],yerr=std_m1[1])
##plt.legend()
#plt.show()
#plt.figure()
#plt.title("distribution of species n2 from "+str(num_traject)+" trajectories with hillfunction g(n2)=q(n2)")
##plt.plot(dist_n[1],label="distribution of species n from "+str(num_traject)+" trajectories")
#plt.errorbar(dist_n2[0],dist_n2[1],yerr=std_n2[1])
##plt.legend()
#plt.figure()
##plt.plot(dist_m[1],label="distribution of species m from "+str(num_traject)+" trajectories")
#plt.title("distribution of species m2 from "+str(num_traject)+" trajectories with hillfunction g(n2)=q(n2)")
#plt.errorbar(dist_m2[0],dist_m2[1],yerr=std_m2[1])
##plt.legend()
#plt.show()

#mean_n1 = 0.
#for i in xrange(0,num_traject):
	#mean_n1 += np.mean(alldata_n1[i].T[1][20000:])
	#hist_i_n1 = np.histogram(alldata_n1[i].T[1][20000:],bins=np.max(alldata_n1[i].T[1][20000:]),range=(0,np.max(alldata_n1[i].T[1][20000:])),normed=True)
	#plt.plot(hist_i_n1[0],label="distribution of species n1 from the "+str(i)+"th trajectory")
#plt.legend()
#mean_n1 = mean_n1/num_traject

#mean_m1 = 0.
#plt.figure()
#for i in xrange(0,num_traject):
	#mean_m1 += np.mean(alldata_m1[i].T[1][20000:])
	#hist_i_m1 = np.histogram(alldata_m1[i].T[1][20000:],bins=np.max(alldata_m1[i].T[1][20000:]),range=(0,np.max(alldata_m1[i].T[1][20000:])),normed=True)
	#plt.plot(hist_i_m1[0],label="distribution of species m1 from the "+str(i)+"th trajectory")
#plt.legend()
#mean_m1 = mean_m1/num_traject

#mean_n2 = 0.
#for i in xrange(0,num_traject):
	#mean_n2 += np.mean(alldata_n2[i].T[1][20000:])
	#hist_i_n2 = np.histogram(alldata_n2[i].T[1][20000:],bins=np.max(alldata_n2[i].T[1][20000:]),range=(0,np.max(alldata_n2[i].T[1][20000:])),normed=True)
	#plt.plot(hist_i_n2[0],label="distribution of species n2 from the "+str(i)+"th trajectory")
#plt.legend()
#mean_n2 = mean_n2/num_traject

#mean_m2 = 0.
#plt.figure()
#for i in xrange(0,num_traject):
	#mean_m2 += np.mean(alldata_m2[i].T[1][20000:])
	#hist_i_m2 = np.histogram(alldata_m2[i].T[1][20000:],bins=np.max(alldata_m2[i].T[1][20000:]),range=(0,np.max(alldata_m2[i].T[1][20000:])),normed=True)
	#plt.plot(hist_i_m2[0],label="distribution of species m from the "+str(i)+"th trajectory")
#plt.legend()
#mean_m2 = mean_m2/num_traject




#def find_histogram(species_num,trajectory,starttime,endtime,smod):
	#'''find histogram of specific trajectory of the stochpy simulation object smod in a certain timeintervall'''
	#smod.GetTrajectoryData(trajectory)
	#data = smod.data_stochsim.getDataInTimeInterval(starttime+0.5*(endtime-starttime),0.5*(endtime-starttime)) # data at each reaction time in an array in form: [[t0,n0,m0],[t1,n1,m1],...]
	#hist = np.histogram(data.T[1],bins=np.max(data.T[species_num]),range=(0,np.max(data.T[species_num])),normed=True)
	#return (hist[1][:-1],hist[0])

#def hist_list(species_num,trajectories,starttime,endtime,smod):
	#''''''
	#alltra = []
	#for tra in np.arange(1,trajectories):
		#hist = find_histogram(species_num,tra,starttime,endtime,smod)[1]
		#alltra.append(hist)
	#maxnum = np.max(map(np.shape,alltra))
	#new_alltra = [np.pad(a, (0,maxnum-np.shape(a)[0]), 'constant') for a in alltra]
	#return new_alltra

#def plot_all_hist(species_num,trajectories,starttime,endtime,smod,text):
	#'''returns distribution of the species specified by its number in the smod.data_stochsim.getDataInTimeInterval(...) array'''
	#new_alltra = hist_list(species_num,trajectories,starttime,endtime,smod)
	#plt.plot(np.sum(new_alltra,axis=0)/np.sum(new_alltra),label=text)
	#return np.sum(new_alltra,axis=0)/np.sum(new_alltra)
	
#def plot_Gill_vs_Spec(trajectories,starttime,endtime,smod,p_n_Spec,p_m_Spec):
	#plt.figure()
	#plt.title('marginal probability density of input proteins n')
	#plt.xlabel('n')
	#plt.ylabel('p(n)')
	#plt.plot(plot_all_hist(1,trajectories,starttime,endtime,smod,'Gillespie'))
	#plt.plot(p_n_Spec,'r',label='Spectral Method')
	#plt.legend()
	#plt.show()

	#plt.figure()
	#plt.title('marginal probability density of output protein m')
	#plt.xlabel('m')
	#plt.ylabel('p(m)')
	#plt.plot(plot_all_hist(2,trajectories,starttime,endtime,smod,'Gillespie'))
	#plt.plot(p_m_Spec,'r',label='Spectral Method')
	#plt.legend()
	#plt.show()

#p_n1 = plot_all_hist(1,10,100,1000,smod,'histogram of species n1')
#np.save('Data/p_n1_10_100_1000.npy',p_n1)

#p_n2 = plot_all_hist(2,10,100,1000,smod,'histogram of species n2')
#np.save('Data/p_n2_10_100_1000.npy',p_n2)

#p_m1 = plot_all_hist(3,10,100,1000,smod,'histogram of species m1')
#np.save('Data/p_m1_10_100_1000.npy',p_m1)

#p_m2 = plot_all_hist(4,10,100,1000,smod,'histogram of species m2')
#np.save('Data/p_m2_10_100_1000.npy',p_m2)


#p_n1_Gill = smod.data_stochsim.species_distributions[0]#shows distribution from last trajectory!
#p_n2_Gill = smod.data_stochsim.species_distributions[1]
#p_m1_Gill = smod.data_stochsim.species_distributions[2]
#p_m2_Gill = smod.data_stochsim.species_distributions[3]
##stochpy.SaveInteractiveSession("interactiveSession1.py")

##plt.figure()
##plt.title('marginal probability density of input proteins n1')
##plt.xlabel('n1')
##plt.ylabel('p(n1)')
##plt.plot(plot_all_hist(1,trajectories,starttime,endtime,smod,'Gillespie'))
##plt.plot(p_n_Spec,'r',label='Spectral Method')
##plt.legend()
##plt.show()

##plt.figure()
##plt.title('marginal probability density of input proteins n1')
##plt.xlabel('n2')
##plt.ylabel('p(n2)')
##plt.plot(plot_all_hist(2,trajectories,starttime,endtime,smod,'Gillespie'))
##plt.plot(p_n_Spec,'r',label='Spectral Method')
##plt.legend()
##plt.show()

##plt.figure()
##plt.title('marginal probability density of output protein m')
##plt.xlabel('m')
##plt.ylabel('p(m)')
##plt.plot(plot_all_hist(3,trajectories,starttime,endtime,smod,'Gillespie'))
##plt.plot(p_m_Spec,'r',label='Spectral Method')
##plt.legend()
##plt.show()

##plt.figure()
##plt.title('marginal probability density of output protein m2')
##plt.xlabel('m2')
##plt.ylabel('p(m2)')
##plt.plot(plot_all_hist(4,trajectories,starttime,endtime,smod,'Gillespie'))
##plt.plot(p_m_Spec,'r',label='Spectral Method')
##plt.legend()
##plt.show()

#plt.figure()
#plt.title('marginal probability density of input proteins n1')
#plt.xlabel('n1')
#plt.ylabel('p(n1)')
#plt.plot(p_n1_Gill[0],p_n1_Gill[1],label='Gillespie')
#plt.plot(p_n1_Spec,'r',label='Spectral Method')
#plt.legend()
#plt.show()

#plt.figure()
#plt.title('marginal probability density of output protein n2')
#plt.xlabel('n2')
#plt.ylabel('p(n2)')
#plt.plot(p_n2_Gill[0],p_n2_Gill[1],label='Gillespie')
#plt.plot(p_n2_Spec,'r',label='Spectral Method')
#plt.legend()
#plt.show()

#plt.figure()
#plt.title('marginal probability density of output protein m1')
#plt.xlabel('m1')
#plt.ylabel('p(m1)')
#plt.plot(p_m1_Gill[0],p_m1_Gill[1],label='Gillespie')
#plt.plot(p_m1_Spec,'r',label='Spectral Method')
#plt.legend()
#plt.show()

#plt.figure()
#plt.title('marginal probability density of output protein m2')
#plt.xlabel('m2')
#plt.ylabel('p(m2)')
#plt.plot(p_m2_Gill[0],p_m2_Gill[1],label='Gillespie')
#plt.plot(p_m2_Spec,'r',label='Spectral Method')
#plt.legend()
#plt.show()