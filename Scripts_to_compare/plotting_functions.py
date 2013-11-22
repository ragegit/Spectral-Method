import numpy as np
import matplotlib.pyplot as plt
import math

# functions to load the full probability distribution found with the Spectral Method from file:
def load_GeneRegulation_Diffusion_SpM(filename):
	p_mat_m2_n2_m1_n1 = np.load('../SpectralMethod/GeneRegulation_Diffusion_SpM/Data/'+filename)
	p_n2_m1_n1 = np.sum(p_mat_m2_n2_m1_n1,axis = 0)
	p_m1_n1 = np.sum(p_n2_m1_n1,axis = 0)
	p_n1 = np.sum(p_m1_n1,axis = 0)
	p_n2 = np.sum(np.sum(p_n2_m1_n1,axis = 1),axis=1)
	p_m1 = np.sum(p_m1_n1,axis = 1)
	p_m2 = np.sum(np.sum(np.sum(p_mat_m2_n2_m1_n1,axis=1),axis=1),axis=1)
	return {'n1':p_n1,'n2':p_n2,'m1':p_m1,'m2':p_m2}# add 'N':N, 'M':M
	
def load_Diffusion_SpM(filename):
	p_mat = np.load('../SpectralMethod/Diffusion_SpM/Data/'+filename)
	p_n1 = np.sum(p_mat, axis=1)
	p_n2 = np.sum(p_mat, axis=0)
	return {'n1':p_n1,'n2':p_n2}

def load_GeneRegulation_SpM(filename):
	p_mat = np.load('../SpectralMethod/GeneRegulation_SpM/Data/'+filename)
	p_n = np.sum(p_mat, axis=1) # sum of all elements in a row n gives value of the marginal distribution p(n)
	p_m = np.sum(p_mat, axis=0)
	return {'n':p_n,'m':p_m}
	
# functions to load the full probability distribution found with Gillespie from file:
def load_GeneRegulation_Diffusion_Gillespie(trajectories,starttime,endtime):
	''' returns the full probability distribution from the two sites Generegulation-Diffusion problem stored in the file named filename as a dictionary. access the marginal distributions via its label: prob_dist = load_GeneRegulation_Diffusion_Gillespie(trajectories,starttime,endtime)'''
	p_n1 = np.load('../Gillespie/GeneRegulation_Diffusion_Gillespie/Data/p_n1_'+str(trajectories)+'_'+str(starttime)+'_'+str(endtime)+'.npy')
	p_n2 = np.load('../Gillespie/GeneRegulation_Diffusion_Gillespie/Data/p_n2_'+str(trajectories)+'_'+str(starttime)+'_'+str(endtime)+'.npy')
	p_m1 = np.load('../Gillespie/GeneRegulation_Diffusion_Gillespie/Data/p_m1_'+str(trajectories)+'_'+str(starttime)+'_'+str(endtime)+'.npy')
	p_m2 = np.load('../Gillespie/GeneRegulation_Diffusion_Gillespie/Data/p_m2_'+str(trajectories)+'_'+str(starttime)+'_'+str(endtime)+'.npy')
	return {'n1':p_n1,'n2':p_n2,'m1':p_m1,'m2':p_m2}

def load_Diffusion_Gillespie(trajectories,starttime,endtime):#??? maybe better to check, if convergence of moments instead of starttime and endtime
	''' returns the full probability distribution from the two sites Generegulation-Diffusion problem stored in the file named filename as a dictionary. access the marginal distributions via its label: prob_dist = load_GeneRegulation_Diffusion_Gillespie(trajectories,starttime,endtime). then prob_dist['n1'] gives the marginal distribution of species n1 '''
	p_n1 = np.load('../Gillespie/Diffusion_Gillespie/Data/p_n1_'+str(trajectories)+'_'+str(starttime)+'_'+str(endtime)+'.npy')
	p_n2 = np.load('../Gillespie/Diffusion_Gillespie/Data/p_n2_'+str(trajectories)+'_'+str(starttime)+'_'+str(endtime)+'.npy')
	return {'n1':p_n1,'n2':p_n2}

def load_GeneRegulation_Gillespie(trajectories,starttime,endtime):
	''' returns the full probability distribution from the two species Generegulation problem stored in the file named filename as a dictionary. access the marginal distributions via its label: prob_dist = load_GeneRegulation_Gillespie(trajectories,starttime,endtime). then prob_dist['n'] gives the marginal distribution of species n'''
	p_n = np.load('../Gillespie/GeneRegulation_Gillespie/Data/p_n_'+str(trajectories)+'_'+str(starttime)+'_'+str(endtime)+'_poisson_g7.npy')
	p_m = np.load('../Gillespie/GeneRegulation_Gillespie/Data/p_m_'+str(trajectories)+'_'+str(starttime)+'_'+str(endtime)+'_poisson_g7.npy')
	return {'n':p_n,'m':p_m}

def plot_Gill_vs_Spec(speciesname,shift1,shift2,p_n_Spec,p_n_Gill):
	plt.figure()
	plt.title('marginal probability density of species '+speciesname+' with shifts g='+str(shift1)+' and q='+str(shift2))
	plt.xlabel(speciesname)
	plt.ylabel('p('+speciesname+')')
	plt.plot(p_n_Gill,'b',label='Gillespie')
	plt.plot(p_n_Spec,'r',label='Spectral Method')
	plt.legend()
	plt.show()

# plots to compare pure gene regulation problem solved with Gillespie and the Spectral Method
#p_n_m_Gillespie = load_GeneRegulation_Gillespie(100,100,1000)
#for j in [40,50]:
	#for g in np.arange(9,12):
		#for q in np.arange(9,12):	
			#p_n_m_SpM = load_GeneRegulation_SpM('p_'+str(j)+'_'+str(j)+'_'+str(j)+'_'+str(j)+'_shifts_'+str(g)+'_'+str(q)+'.npy')
			#plot_Gill_vs_Spec('n',str(g),str(q),p_n_m_SpM['n'],p_n_m_Gillespie['n'])
			#plot_Gill_vs_Spec('m',str(g),str(q),p_n_m_SpM['m'],p_n_m_Gillespie['m'])

p_n_m_Gillespie = load_GeneRegulation_Gillespie(100,100,1000)
for j in [40,50,60]:
	p_n_m_SpM = load_GeneRegulation_SpM('initial_hill3/p_'+str(j)+'_'+str(j)+'_'+str(j)+'_'+str(j)+'_shifts_average.npy')
	f, axarr = plt.subplots(2,1)
	axarr[0].plot(p_n_m_Gillespie['n'],label='species n')
	axarr[1].plot(p_n_m_Gillespie['m'],label='species m')
	axarr[0].plot(p_n_m_SpM['n'],label='species n')
	axarr[0].set_title('species n')
	axarr[1].plot(p_n_m_SpM['m'],label='species m')
	axarr[1].set_title('species m')
	plt.savefig('Plots/GeneRegulation/initial_poisson_g7/p_four_subplots_Genereg_SpM_vs_Gillespie_N'+str(j)+'.png')


p_n_m_Gillespie = load_GeneRegulation_Gillespie(100,100,1000)
for j in [40,50,60]:
	for g in np.arange(9,12):
		for q in np.arange(9,12):	
			p_n_m_SpM = load_GeneRegulation_SpM('initial_hill2/p_'+str(j)+'_'+str(j)+'_'+str(j)+'_'+str(j)+'_shifts_'+str(g)+'_'+str(q)+'.npy')
			f, axarr = plt.subplots(2,1)
			axarr[0].plot(p_n_m_Gillespie['n'],label='species n, g='+str(g)+' q='+str(q))
			axarr[1].plot(p_n_m_Gillespie['m'],label='species m, g='+str(g)+' q='+str(q))
			axarr[0].plot(p_n_m_SpM['n'],label='species n, g='+str(g)+' q='+str(q))
			axarr[0].set_title('species n, g='+str(g)+' q='+str(q))
			axarr[1].plot(p_n_m_SpM['m'],label='species m, g='+str(g)+' q='+str(q))
			axarr[1].set_title('species m, g='+str(g)+' q='+str(q))
			plt.savefig('Plots/GeneRegulation/p_four_subplots_Genereg_SpM_vs_Gillespie_N'+str(j)+'_g'+str(g)+'_q'+str(q)+'.png')


# plots to compare gene regulation and diffusion problem solved with Gillespie and the Spectral Method
#p_n1_n2_m1_m2_Gillespie = load_GeneRegulation_Diffusion_Gillespie(10,100,1000)
#for j in [20,25]:
	#for g in np.arange(9,10):
		#for q in np.arange(9,13):	
			#p_n1_n2_m1_m2_SpM = load_GeneRegulation_Diffusion_SpM('initial1/p_'+str(j)+'_'+str(j)+'_'+str(j)+'_'+str(j)+'_shifts_'+str(g)+'_'+str(q)+'.npy')
			#plot_Gill_vs_Spec('n1',str(g),str(q),p_n1_n2_m1_m2_SpM['n1'],p_n1_n2_m1_m2_Gillespie['n1'])
			#plot_Gill_vs_Spec('m1',str(g),str(q),p_n1_n2_m1_m2_SpM['m1'],p_n1_n2_m1_m2_Gillespie['m1'])
			#plot_Gill_vs_Spec('n2',str(g),str(q),p_n1_n2_m1_m2_SpM['n2'],p_n1_n2_m1_m2_Gillespie['n2'])
			#plot_Gill_vs_Spec('m2',str(g),str(q),p_n1_n2_m1_m2_SpM['m2'],p_n1_n2_m1_m2_Gillespie['m2'])

p_n1_n2_m1_m2_Gillespie = load_GeneRegulation_Diffusion_Gillespie(10,100,1000)
for j in [15,20,25]:
	for g in np.arange(9,13):
		for q in np.arange(9,13):	
			p_n1_n2_m1_m2_SpM = load_GeneRegulation_Diffusion_SpM('initial4/p_'+str(j)+'_'+str(j)+'_'+str(j)+'_'+str(j)+'_shifts_'+str(g)+'_'+str(q)+'.npy')
			f, axarr = plt.subplots(2, 2)
			axarr[0,0].plot(p_n1_n2_m1_m2_Gillespie['n1'],label='species n1, g='+str(g)+' q='+str(q))
			axarr[0,1].plot(p_n1_n2_m1_m2_Gillespie['m1'],label='species m1, g='+str(g)+' q='+str(q))
			axarr[1,0].plot(p_n1_n2_m1_m2_Gillespie['n2'],label='species n2, g='+str(g)+' q='+str(q))
			axarr[1,1].plot(p_n1_n2_m1_m2_Gillespie['m2'],label='species m2, g='+str(g)+' q='+str(q))
			
			axarr[0,0].plot(p_n1_n2_m1_m2_SpM['n1'],label='species n1, g='+str(g)+' q='+str(q))
			axarr[0,0].set_title('species n1, g='+str(g)+' q='+str(q))
			axarr[0,1].plot(p_n1_n2_m1_m2_SpM['m1'],label='species m1, g='+str(g)+' q='+str(q))
			axarr[0,1].set_title('species m1, g='+str(g)+' q='+str(q))
			axarr[1,0].plot(p_n1_n2_m1_m2_SpM['n2'],label='species n2, g='+str(g)+' q='+str(q))
			axarr[1,0].set_title('species n2, g='+str(g)+' q='+str(q))
			axarr[1,1].plot(p_n1_n2_m1_m2_SpM['m2'],label='species m2, g='+str(g)+' q='+str(q))
			axarr[1,1].set_title('species m2, g='+str(g)+' q='+str(q))
			plt.savefig('Plots/GeneRegulation_Diffusion/initial4/p_four_subplots_Genereg_Diff_SpM_vs_Gillespie_N'+str(j)+'_g'+str(g)+'_q'+str(q)+'.png')

p_n1_n2_m1_m2_Gillespie = load_GeneRegulation_Diffusion_Gillespie(10,100,1000)
for j in [20,25]:
	for g in np.arange(9,13):
		for q in np.arange(9,13):	
			p_n1_n2_m1_m2_SpM = load_GeneRegulation_Diffusion_SpM('initial2/p_'+str(j)+'_'+str(j)+'_'+str(j)+'_'+str(j)+'_shifts_'+str(g)+'_'+str(q)+'.npy')
			f, axarr = plt.subplots(2, 2)
			axarr[0,0].plot(p_n1_n2_m1_m2_Gillespie['n1'],label='species n1, g='+str(g)+' q='+str(q))
			axarr[0,1].plot(p_n1_n2_m1_m2_Gillespie['m1'],label='species m1 g='+str(g)+' q='+str(q))
			axarr[1,0].plot(p_n1_n2_m1_m2_Gillespie['n2'],label='species n2 g='+str(g)+' q='+str(q))
			axarr[1,1].plot(p_n1_n2_m1_m2_Gillespie['m2'],label='species m2 g='+str(g)+' q='+str(q))
			
			axarr[0,0].plot(p_n1_n2_m1_m2_SpM['n1'],label='species n1, g='+str(g)+' q='+str(q))
			axarr[0,0].set_title('species n1, g='+str(g)+' q='+str(q))
			axarr[0,1].plot(p_n1_n2_m1_m2_SpM['m1'],label='species m1, g='+str(g)+' q='+str(q))
			axarr[0,1].set_title('species m1, g='+str(g)+' q='+str(q))
			axarr[1,0].plot(p_n1_n2_m1_m2_SpM['n2'],label='species n2, g='+str(g)+' q='+str(q))
			axarr[1,0].set_title('species n2, g='+str(g)+' q='+str(q))
			axarr[1,1].plot(p_n1_n2_m1_m2_SpM['m2'],label='species m2, g='+str(g)+' q='+str(q))
			axarr[1,1].set_title('species m2, g='+str(g)+' q='+str(q))
			plt.savefig('Plots/GeneRegulation_Diffusion/initial2/p_four_subplots_Genereg_Diff_SpM_vs_Gillespie_N'+str(j)+'_g'+str(g)+'_q'+str(q)+'.png')

#f, axarr = plt.subplots(3, 1)
#for g in np.arange(9,12):
	##for q in np.arange(9,12):
	#p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(9)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	#p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	#p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	#p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	#p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	#p_m1 = np.sum(p_m2_m1, axis = 1)
	#p_m2 = np.sum(p_m2_m1, axis = 0)
	#axarr[0].plot(p_m1,label='g='+str(g) )
	#axarr[0].set_title('q=9')
	#p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(10)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	#p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	#p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	#p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	#p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	#p_m1 = np.sum(p_m2_m1, axis = 1)
	#p_m2 = np.sum(p_m2_m1, axis = 0)
	#axarr[1].plot(p_m1,label='g='+str(g) )
	#axarr[1].set_title('q=10')
	#p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(11)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	#p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	#p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	#p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	#p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	#p_m1 = np.sum(p_m2_m1, axis = 1)
	#p_m2 = np.sum(p_m2_m1, axis = 0)
	#axarr[2].plot(p_m1,label='g='+str(g) )
	#axarr[2].set_title('q=11')
#plt.legend()
#plt.ylabel('marginal probability density p(m1)')
#plt.xlabel('protein number m1')
#plt.show()