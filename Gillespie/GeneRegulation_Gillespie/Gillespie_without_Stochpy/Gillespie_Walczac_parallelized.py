import numpy as np
from initializations_Gillespie import*
from functions import *
import matplotlib.pyplot as plt

def Gillespie_algo(t_max,nu):
	num_traj = 2**nu
	num_t_save = np.int(num_traj/1000)
	tract_storage = np.zeros((num_traj,num_t_save,2)) # 3D array to store "num_traj" numbers of trajectories where for each trajectory every 1000^th timestep the number of particles of each species is saved
	tract_storage[:,0,0] = 10 # initial number of particles of type n. Could also be initialized with a random number
	tract_storage[:,0,1] = 0 # initial number of particles of type m
	rates = np.zeros((num_traj,4)) # initial rates for each trajectory
	rates[:,0] = tract_storage[:,0,0]
	rates[:,1] = tract_storage[:,0,1]
	rates[:,2] = g
	rates[:,3] = map(q_n0,tract_storage[:,0,1])
	for t in np.arange(1,t_max):
		r_tau = 1. - np.random.rand(1) # random number to find time at which the next reaction occurs
		react = np.random.rand(1) # random number to find which reaction occurs next
		rate_sum_total = np.sum(rates) # the sum of all rates of all trajectories as used in a general Gillespie algorithm
		
		# find out for which trajectory the reaction occurs:
		new_rate_sum = np.sum(rates, axis=1)
		new_interval_length = np.sum(new_rate_sum[:(2**(nu-1))])
		mini = 2**0
		maxi = 2**nu
		for i in np.arange(1,nu):
			if( react*new_interval_length < new_interval_length):
				new_rate_sum = new_rate_sum[:(2**(nu-i))]
				new_interval_length = np.sum(new_rate_sum)
				maxi = maxi - 2**(nu-i)
			else:
				new_rate_sum = new_rate_sum[(2**(nu-i)):]
				new_interval_length = np.sum(new_rate_sum)
				counter = counter*3/4
				mini = mini + 2**(nu-i)
			








def Gillespie_algo(reacmax):
	'''function which draws a trajectory with "reacmax" reactions using a Gillespie algorithm for the stoch. diff. eq.:
	p(n,m) = g*p(n-1,m) + (n+1)*p(n+1,m) - (g+n)*p(n,m) + rho*(q(n,n0)*p(n,m-1) + (m+1)*p(n,m+1) - (q(n,n0)+m)*p(n,m)) '''
	#initializations for each trajectory:
	rho = 1. # or is rho different?
	n_ini = np.random.randint(0,nmax) # initial number of particles of type n
	m_ini = np.random.randint(0,mmax) # initial number of particles of type m
	deg_n = 1.*n_ini # initial degradation rate of species n
	deg_m = rho*m_ini # initial degradation rate of species m
	cre_n = g # constant creation rate of species n
	cre_m = q(n_ini,n0) # n dependend creation rate of species m
	rates = np.array([deg_n, deg_m, cre_n, cre_m]) # deg_n = n+g , deg_m = (m + q(n,n0)))*rho , cre_n = (n-1)+g , cre_m = ((m-1)+q(n,n0))*rho
	particle_storage = np.zeros((3,reacmax))
	particle_storage[0][0] = 0.
	particle_storage[1][0] = n_ini
	particle_storage[2][0] = m_ini
	for t in np.arange(1,reacmax):
		r_tau = 1. - np.random.rand(1) # random number to find time at which the next reaction occurs
		react = np.random.rand(1) # random number to find which reaction occurs next
		rate_sum = np.float(np.sum(rates)) # the sum of all rates as used in a general Gillespie algo
		#if(rate_sum == 0):
			#print(rates)
			#print(t)
			#print(particle_storage)
			#print(rate_sum)
		tau = -1./rate_sum*np.log(r_tau) # time until next reaction occurs
		# updates:
		# total time update:
		particle_storage[0][t] = particle_storage[0][t-1]+tau
		# particle number and rate updates:
		if (react < np.float(rates[0])/rate_sum and react > 0.):
			particle_storage[1][t] = particle_storage[1][t-1]-1.
			rates[0] = particle_storage[1][t]
			rates[3] = q(particle_storage[1][t],n0) # n dependend creation rate of species m
		elif (react < np.float(rates[0]+rates[1])/rate_sum and react > np.float(rates[0])/rate_sum):
			particle_storage[2][t] = particle_storage[2][t-1]-1.
			rates[1] = rates[1]-rho
		elif (react < np.float(rates[0]+rates[1]+rates[2])/rate_sum and react > np.float(rates[0]+rates[1])/rate_sum):
			particle_storage[1][t] = particle_storage[1][t-1]+1.
			rates[0] = rates[0] + 1.
			rates[3] = q(particle_storage[1][t],n0) # n dependend creation rate of species m
		elif (react < 1. and react > np.float(rates[0]+rates[1]+rates[2])/rate_sum):
			particle_storage[2][t] = particle_storage[2][t-1]+1
			rates[1] = rates[1]+rho
	return particle_storage