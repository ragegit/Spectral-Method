import numpy as np
from initializations_Gillespie import *
from rate_functions import *
import matplotlib.pyplot as plt

def Gillespie_algo(reacmax):
	'''function which draws a trajectory with "reacmax" reactions using a Gillespie algorithm for the stoch. diff. eq.:
	p(n,m) = g*p(n-1,m) + (n+1)*p(n+1,m) - (g+n)*p(n,m) + rho*(q(n,n0)*p(n,m-1) + (m+1)*p(n,m+1) - (q(n,n0)+m)*p(n,m)) '''
	
	# initializations for each trajectory:
	rho = 1. # or is rho different?
	n_ini = 1. #np.random.randint(0,nmax) # initial number of particles of type n
	m_ini = 1. #np.random.randint(0,mmax) # initial number of particles of type m
	deg_n = 1.*n_ini # initial degradation rate of species n
	deg_m = rho*m_ini # initial degradation rate of species m
	cre_n = g_n0(n_ini) # n dependend creation rate of species n (autoregulation)
	cre_m = q_n0(n_ini) # n dependend creation rate of species m (regulated by n)
	rates = np.array([deg_n, deg_m, cre_n, cre_m]) # deg_n = n+g , deg_m = (m + q(n,n0)))*rho , cre_n = (n-1)+g , cre_m = ((m-1)+q(n,n0))*rho
	particle_storage = np.zeros((3,reacmax)) # array, which stores the time a reaction occurs and the number of particles at each such time
	particle_storage[0,0] = 0.
	particle_storage[1,0] = n_ini
	particle_storage[2,0] = m_ini
	
	for t in np.arange(reacmax-1):
		# calculate end of each rate interval
		rate_sum = np.sum(rates) # the sum of all rates as used in a general Gillespie algo
		rate_sum_deg_n = np.sum(rates[:1]) # end of the first rate interval
		rate_sum_deg_m = np.sum(rates[:2]) # end of the second rate interval
		rate_sum_cre_n = np.sum(rates[:3]) # end of the third rate interval
		
		r_tau = 1. - np.random.rand(1) # random number to find time at which the next reaction occurs
		react = np.random.rand(1) # random number to find which reaction occurs next
		tau = -1./rate_sum*np.log(r_tau) # time until next reaction occurs
		# updates:
		# total time update:
		particle_storage[0,t+1] = particle_storage[0,t]+tau
		#print((particle_storage[0,t],particle_storage[1,t],particle_storage[2,t]))
		#print(rates)
		# particle number and rate updates:
		if (rate_sum*react <= rate_sum_deg_n and rate_sum*react >= 0.):
			particle_storage[1,t+1] = particle_storage[1,t]-1.
			particle_storage[2,t+1] = particle_storage[2,t]
			rates[0] = particle_storage[1,t+1]
			rates[2] = g_n0(particle_storage[1,t+1]) # n dependend creation rate of species n
			rates[3] = q_n0(particle_storage[1,t+1]) # n dependend creation rate of species m
		elif (rate_sum*react <= rate_sum_deg_m and rate_sum*react > rate_sum_deg_n):
			particle_storage[2,t+1] = particle_storage[2,t]-1.
			particle_storage[1,t+1] = particle_storage[1,t]
			rates[1] = particle_storage[2,t+1]*rho
		elif (rate_sum*react <= rate_sum_cre_n and rate_sum*react > rate_sum_deg_m): 
			particle_storage[1,t+1] = particle_storage[1,t]+1.
			particle_storage[2,t+1] = particle_storage[2,t]
			rates[0] = particle_storage[1,t+1]
			rates[2] = g_n0(particle_storage[1,t+1]) # n dependend creation rate of species n
			rates[3] = q_n0(particle_storage[1,t+1]) # n dependend creation rate of species m
		elif (rate_sum*react <= rate_sum and rate_sum*react > rate_sum_cre_n):
			particle_storage[2,t+1] = particle_storage[2,t]+1
			particle_storage[1,t+1] = particle_storage[1,t]
			rates[1] = particle_storage[2,t+1]
	return particle_storage