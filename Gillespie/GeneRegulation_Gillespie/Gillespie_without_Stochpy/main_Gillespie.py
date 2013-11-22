from Gillespie_Walczac import *

''' sample from trajectories after the steady state has been reached
two possibilities of sampling due to statistical independence:
- sample many times from one trajectory as soon as the system is in steady state
- sample once for each trajectory to reasure results are independent '''

reacmax = 10000.
max_trajectory = 10
tract_storage = [] #np.zeros(max_trajectory)

for i in np.arange(max_trajectory):
	#np.random.seed(i)
	next_trajectory = Gillespie_algo(reacmax)
	tract_storage.append(next_trajectory)

# plot single trajectory
plt.figure()
plt.xlabel('time')
plt.ylabel('number of proteins of input species n')
plt.plot(tract_storage[9][0,:],tract_storage[9][1,:])
plt.figure()
plt.xlabel('time')
plt.ylabel('number of proteins of output species m')
plt.plot(tract_storage[9][0,:],tract_storage[9][2,:])
plt.show()

#plot histogram of one trajectory
plt.figure()
plt.title('histogram for a single trajectory')
plt.xlabel('number of proteins of input species n')
plt.ylabel('histogram for the marginal distribution p(n)')
plt.hist(tract_storage[9][1,1000:], bins=nmax, range=(0,nmax))
plt.savefig('histogram_single_traj_n.png')
plt.figure()
plt.xlabel('number of proteins of output species m')
plt.ylabel('histogram for the marginal distribution p(m)')
plt.hist(tract_storage[9][2,1000:], bins=mmax, range=(0,mmax))
plt.savefig('histogram_single_traj_m.png')

# draw histogram using the data from all trajectories
all_traj_equi_n = np.zeros((reacmax-reacmax/max_trajectory)*max_trajectory)
all_traj_equi_m = np.zeros((reacmax-reacmax/max_trajectory)*max_trajectory)
j=0
for element in tract_storage:
	k=j+(max_trajectory-1)
	all_traj_equi_n[(j*reacmax/max_trajectory):(k*reacmax/max_trajectory)] = element[1][reacmax/max_trajectory:]
	all_traj_equi_m[(j*reacmax/max_trajectory):(k*reacmax/max_trajectory)] = element[2][reacmax/max_trajectory:]
	j=k


#plot histogram for data from all trajectories
plt.figure()
plt.title('histogram found from many reactions at steady state and many trajectories')
plt.xlabel('number of proteins of input species n')
plt.ylabel('histogram for the marginal distribution p(n)')
plt.hist(all_traj_equi_n, bins=nmax, range=(0,nmax))
plt.savefig('histogram_1000_traj_n.png')
plt.figure()
plt.xlabel('number of proteins of output species m')
plt.ylabel('histogram for the marginal distribution p(m)')
plt.hist(all_traj_equi_m, bins=mmax, range=(0,mmax))
plt.savefig('histogram_1000_traj_m.png')

# mean of all trajectories assuming the equilibrium is reached after 1000 reactions
mean_list_n = []
mean_list_m = []

# check if equilibrium is reached after 1000 reactions
for i in range((np.int(reacmax/1000) - 1)):
	equi = i*1000
	collect_n = np.zeros(max_trajectory)
	collect_m = np.zeros(max_trajectory)
	j=0
	for element in tract_storage:
		collect_n[j] = np.mean(element[1][equi:])
		collect_m[j] = np.mean(element[2][equi:])
		j+=1
		
	mean_list_n.append(np.mean(collect_n)) # dictionary would be nicer!!!
	mean_list_m.append(np.mean(collect_m))


mean_list_n = []
mean_list_m = []

for element in tract_storage:
	mean_list_n.append(np.mean(element[1][100:]))
	mean_list_m.append(np.mean(element[2][100:]))
mean_n_100 = np.mean(mean_list_n)
mean_m_100 = np.mean(mean_list_m)

#steady_n = []
#steady_m = []
mean_list_n = []
mean_list_m = []

for element in tract_storage:
	mean_list_n.append(np.mean(element[1][50:]))
	mean_list_m.append(np.mean(element[2][50:]))
mean_n_50 = np.mean(mean_list_n)
mean_m_50 = np.mean(mean_list_m)
