Data is saved as numpy arrays which can be loaded e.g. with a relative path from the folder 'SpectralMethod_vs_Gillespie' via

Some distributions found from the original trajectories is saved in the following format:
p_n = np.save('/Gillespie/GeneRegulation_Gillespie/Data/p_n_'+str(trajectories)+str(starttime)+str(endtime)+'.npy',p_n)

Distributions and standarddeviations from stochpy are saved in the format:
p_n = np.save('/Gillespie/GeneRegulation_Gillespie/Data/p_n_'+str(trajectories)+str(endtime)+'.npy',dist_n)
std_n = np.save('/Gillespie/GeneRegulation_Gillespie/Data/p_n_'+str(trajectories)+str(endtime)+'.npy',std_n)

where trajectories is the number of trajectories, starttime is the time after which data was expected to be at steady state and was counted in the histogram and endtime is the endtime of the simulation
