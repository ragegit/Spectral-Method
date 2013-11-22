import numpy as np
import math
from initializations import *
from functions import *
import numpy.linalg as la
import matplotlib.pyplot as plt

def main_g_q_average(nmax,mmax,jmax,kmax):
	# vector of initial marginal distribution p(n)
	p_vec = p_vector(nmax)

	# vectors of creation rates g(n) and q(n) to find the mean rates for the shifts -> why needed?
	q_vec = rate_vec(nmax,q_n0)
	g_vec = rate_vec(nmax,g_n0)

	q_mean = np.dot(q_vec,p_vec) # constant part of the creation rate of species m
	g_mean = np.dot(g_vec,p_vec) # special in this example: creation rate of both species are regulated by species n with the same functional dependency
	delta_q = delta_matrix(q_n0,nmax,q_mean) # diagonal matrix with entries (q_mean - q(n)) on the n^th diagonal element
	delta_g = delta_matrix(g_n0,nmax,g_mean)
	dual_base_mat_n = dual_mat(g_mean,nmax,jmax) # contains elements <j|n> TODO: not normalized! 
	base_mat_n = mat(g_mean,nmax,jmax) # contains elements <n|j> - normalized :-)
	base_mat_m = mat(q_mean,mmax,kmax) # contains elements <m|k> - normalized :-) 
	new_delta_q = np.dot(np.dot(dual_base_mat_n,delta_q),base_mat_n) # delta_q in new base {|j>}
	new_delta_g = np.dot(np.dot(dual_base_mat_n,delta_g),base_mat_n)
	sub = subdiag(jmax,jmax)

	# gen will contain the coefficients of the generating function in the |j,m> base
	gen = np.zeros((jmax,kmax))

	# first generating function: vector over indecies j with k=0
	gen[:,0] = np.dot(dual_base_mat_n, p_vec)
	allcond = []
	
	for k in np.arange(1,kmax):
		diag = np.zeros((jmax,jmax))
		for j in np.arange(jmax):
			diag[j,j] = j + rho*k
		A = -rho*(diag + np.dot(sub,new_delta_g))
		gen[:,k] = la.lstsq(A, np.dot(new_delta_q,gen[:,k-1]))[0]#gen[:,k] = la.solve(A, np.dot(new_delta_q,gen[:,k-1]))
		cond = la.norm(la.inv(A))*la.norm(A)
		allcond.append(cond)
		#print( np.all(np.dot(new_delta_q,gen[:,k-1]) == np.dot(A,gen[:,k]) ) )
	p_mat = np.dot(base_mat_n,np.dot(gen,base_mat_m.T)) # matrix containing the full distribution, i.e. for row n and column m it contains p(n,m)
	return (p_mat,allcond,g_mean,q_mean)

N = 50
p_mat =  main_g_q_average(N,N,N,N)

p_n = np.sum(p_mat[0], axis=1) # sum of all elements in a row n gives value of the marginal distribution p(n)
p_m = np.sum(p_mat[0], axis=0) # sum of all elements in a column m gives value of the marginal distribution p(m)

n_mean = np.dot(np.arange(N),p_n) # compare with theoretically expected value: <n> = <g(n)>/1 where 1 is the constant degradation rate here
m_mean = np.dot(np.arange(N),p_m)

p_n_m_Gill = np.load('../../Gillespie/GeneRegulation_Gillespie/Data/initial_poisson_g7/p_n_m_1trajectory_endtime500000.npy')
p_n_m_Gill = p_n_m_Gill/np.sum(p_n_m_Gill)
p_n_Gill = np.sum(p_n_m_Gill,axis=1)
p_m_Gill = np.sum(p_n_m_Gill,axis=0)


for j in [N]:
	f, axarr = plt.subplots(2,sharex=True)
	axarr[0].plot(p_m,label='g='+str(p_mat[2])+',q='+str(p_mat[3]))
	axarr[0].plot(p_m_Gill,label='Gillespie')
	axarr[0].set_title('species m for mmax=50 with shifts g='+str(p_mat[2])+'and q='+str(p_mat[3]))
	axarr[1].plot(p_n,label='g='+str(p_mat[2])+'and q='+str(p_mat[3]))
	axarr[1].plot(p_n_Gill,label='Gillespie')
	axarr[1].set_title('species n for nmax=50 with shifts g='+str(p_mat[2])+'and q='+str(p_mat[3]))
	axarr[0].legend()
	axarr[1].legend()
	plt.show()


#for j in [20,30,40,50,60]:
			#p_mat = main_g_q_average(j,j,j,j)
			#np.save('Data/initial_hill3/p_'+str(j)+'_'+str(j)+'_'+str(j)+'_'+str(j)+'_shifts_average.npy',p_mat)

#examine condition of matrices:
def spec_meth(shift_g,shift_q):
	p_vec = p_vector(nmax)
	q_mean = shift_q # constant part of the creation rate of species m
	g_mean = shift_g # special in this example: creation rate of both species are regulated by species n with the same functional dependency
	delta_q = delta_matrix(q_n0,nmax,q_mean) # diagonal matrix with entries (q_mean - q(n)) on the n^th diagonal element
	delta_g = delta_matrix(g_n0,nmax,g_mean)
	dual_base_mat_n = dual_mat(g_mean,nmax,jmax) # contains elements <j|n> TODO: not normalized! 
	base_mat_n = mat(g_mean,nmax,jmax) # contains elements <n|j> - normalized :-)
	base_mat_m = mat(q_mean,mmax,kmax) # contains elements <m|k> - normalized :-) 
	new_delta_q = np.dot(np.dot(dual_base_mat_n,delta_q),base_mat_n) # delta_q in new base {|j>}
	new_delta_g = np.dot(np.dot(dual_base_mat_n,delta_g),base_mat_n)
	sub = subdiag(jmax,jmax)

	# gen will contain the coefficients of the generating function in the |j,m> base
	gen = np.zeros((jmax,kmax))

	# first generating function: vector over indecies j with k=0
	gen[:,0] = np.dot(dual_base_mat_n, p_vec)
	allcond = []
	maxerr = []
	for k in np.arange(1,kmax):
		diag = np.zeros((jmax,jmax))
		for j in np.arange(jmax):
			diag[j,j] = j + rho*k
		#A = -rho*(diag + np.dot(sub,new_delta_g))
		#gen[:,k] = np.dot(la.inv(A), np.dot(new_delta_q,gen[:,k-1])) #la.lstsq(A, np.dot(new_delta_q,gen[:,k-1]))[0]
		A = -rho*(diag + np.dot(sub,new_delta_g))
		gen[:,k] = la.lstsq(A, np.dot(new_delta_q,gen[:,k-1]))[0]#gen[:,k] = la.solve(A, np.dot(new_delta_q,gen[:,k-1]))
		cond = la.norm(la.inv(A))*la.norm(A)
		allcond.append(cond)
		err = np.max(np.abs(np.dot(new_delta_q,gen[:,k-1]) - np.dot(A,gen[:,k])))
		if err>10**(-7):
			#print((err, shift_g, shift_q))
			maxerr.append(err)
		#print( np.all(np.dot(new_delta_q,gen[:,k-1]) == np.dot(A,gen[:,k]) ) )
	p_mat = np.dot(base_mat_n,np.dot(gen,base_mat_m.T)) # matrix containing the full distribution, i.e. for row n and column m it contains p(n,m)
	return [p_mat,allcond,maxerr,shift_g,shift_q,gen]

cond_shift_table = ['(conditions, error, shift g, shift q)']
for g in [7.5,8.,8.5,9.,9.5,9.85,10.,10.5,11.,11.5,12.]:
	for q in [7.5,8.,8.5,9.,9.5,9.85,10.,10.5,11.,11.5,12.]:
		tup = spec_meth(g,q)
		if tup[2]==[]:
			cond_shift_table.append([max(tup[1]),0,tup[3],tup[4]])
		else:
			cond_shift_table.append([max(tup[1]),max(tup[2]),tup[3],tup[4]])
#np.save('Data/conditions_vs_shifts',cond_shift_table)
np.savetxt('Data/conditions_vs_shifts.csv',np.array(cond_shift_table[1:]),fmt='%.4e', delimiter=' ', newline='\n', header='(conditions, error, shift g, shift q)', footer='', comments='# ')

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 9.5}

plt.rc('font', **font)

# Examine the generating function. The procedure is the following, for different shifts, i.e. different base and dual base transformations, find the generating function from the Gillespie simulated p(n,m) and compare it with the generating function found by the spectral method recursion steps.

# load Gillespie p(n,m) from file for poissonian input with g=7 and q(h) a hill function.
p_n_m_Gill = np.load('../../Gillespie/GeneRegulation_Gillespie/Data/initial_poisson_g7/p_n_m_1trajectory_endtime100000.npy')
p_n_m_Gill = p_n_m_Gill/np.sum(p_n_m_Gill)

# enlargen probability function with nullspace
p_n_m_Gill_new = np.zeros((50,50))
p_n_m_Gill_new[:22,:22] = p_n_m_Gill

# compare generating function from Gillespie with that from the spectral method by plotting
f, axarr = plt.subplots(2,sharex=False)#True) #
colorcount = {0:'b',1:'g',2:'r',3:'c',4:'m',5:'y',6:'k',7:'w'}
for j in [50]:
	allcondmax = []
	all_p_mat = []
	i = 0
	for g in [6.0,9.95]:#np.arange(7,12):
		for q in [6,9.95]:#np.arange(7,12):
			(p_mat_n_m,cond,gen) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q),spec_meth(g,q)[-1])#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			dual_Gill_n = dual_mat(g,50,50)
			dual_Gill_m = dual_mat(q,50,50)
			new_gen = np.dot(dual_Gill_n,np.dot(p_mat_n_m[0],dual_Gill_m.T))
			gen_Gill = np.dot(dual_Gill_n,np.dot(p_n_m_Gill_new,dual_Gill_m.T))
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(gen,axis=0),label='g='+str(g)+',q='+str(q),colorcount[i])
			axarr[0].plot(np.sum(gen_Gill,axis=0),'o',label='Gillespie',colorcount[i])
			axarr[0].set_title('species m for mmax=50')
			axarr[1].plot(np.sum(gen,axis=1),label='g='+str(g)+',q='+str(q),colorcount[i])
			axarr[1].plot(np.sum(gen_Gill,axis=1),'o',label='Gillespie',colorcount[i])
			axarr[1].set_title('species n for nmax=50')
			i = i+1
	axarr[0].legend()
	axarr[1].legend()
	plt.xlim(xmax=50)
	plt.show()

#examine the condition of inverted matrices
p_n_m_Gill = np.load('../../Gillespie/GeneRegulation_Gillespie/Data/initial_hill2/p_n_m_1trajectory_endtime500000.npy')
p_n_m_Gill = p_n_m_Gill/np.sum(p_n_m_Gill)

p_m_Gill = np.sum(p_n_m_Gill,axis=0)
p_n_Gill = np.sum(p_n_m_Gill,axis=1)
# shifts0
f, axarr = plt.subplots(2,sharex=True)
#colorcount = {0:'b',1:'g',2:'r',3:'c',4:'m',5:'y',6:'k',7:'w'}
for j in [50]:
	allcondmax = []
	all_p_mat = []
	for g in [6.0,9.86,10.5]:#np.arange(7,12):
		for q in [8.0,9.86,9.9]:#np.arange(7,12):
			(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
			#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
			axarr[0].set_title('p(m) as function of particle number m for all cutoffs 50')
			axarr[1].plot(np.sum(p_mat_n_m[0],axis=1),label='g='+str(g)+',q='+str(q))
			axarr[1].set_title('p(n) as function of particle number n for all cutoffs 50')
	axarr[0].plot(p_m_Gill,'o',label='Gillespie')
	axarr[0].legend()
	axarr[1].plot(p_n_Gill,'o',label='Gillespie')
	axarr[1].legend()
	plt.xlim(xmax=30)
	plt.show()

# shifts1
f, axarr = plt.subplots(2,sharex=True)
for j in [50]:
	allcondmax = []
	all_p_mat = []
	for g in [7.5,8.,8.5,9.]:#np.arange(7,12):
		for q in [7.5,8.,8.5,9.]:#np.arange(7,12):
			(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
			#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
			axarr[0].set_title('species m for mmax=50')
			axarr[1].plot(np.sum(p_mat_n_m[0],axis=1),label='g='+str(g)+',q='+str(q))
			axarr[1].set_title('species n for nmax=50')
	axarr[0].legend()
	axarr[1].legend()
	plt.xlim(xmax=30)
	plt.show()

# shifts2
f, axarr = plt.subplots(2,sharex=True)
for j in [50]:
	allcondmax = []
	all_p_mat = []
	for g in [9.5,10.,10.5,11.]:#np.arange(7,12):
		for q in [7.5,8.,8.5,9.]:#np.arange(7,12):
			(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
			#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
			axarr[0].set_title('species m for mmax=50')
			axarr[1].plot(np.sum(p_mat_n_m[0],axis=1),label='g='+str(g)+',q='+str(q))
			axarr[1].set_title('species n for nmax=50')
	axarr[0].legend()
	axarr[1].legend()
	plt.xlim(xmax=30)
	plt.show()


# shifts3
f, axarr = plt.subplots(2,sharex=True)
for j in [50]:
	allcondmax = []
	all_p_mat = []
	for g in [7.5,8.,8.5,9.]:#np.arange(7,12):
		for q in [9.5,10.,10.5,11.]:#np.arange(7,12):
			(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
			#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
			axarr[0].set_title('species m for mmax=50')
			axarr[1].plot(np.sum(p_mat_n_m[0],axis=1),label='g='+str(g)+',q='+str(q))
			axarr[1].set_title('species n for nmax=50')
	axarr[0].legend()
	axarr[1].legend()
	plt.xlim(xmax=30)
	plt.show()

# shifts4
f, axarr = plt.subplots(2,sharex=True)
for j in [50]:
	allcondmax = []
	all_p_mat = []
	for g in [9.5,10.,10.5,11.]:#np.arange(7,12):
		for q in [9.5,10.,10.5,11.]:#np.arange(7,12):
			(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
			#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
			axarr[0].set_title('species m for mmax=50')
			axarr[1].plot(np.sum(p_mat_n_m[0],axis=1),label='g='+str(g)+',q='+str(q))
			axarr[1].set_title('species n for nmax=50')
	axarr[0].legend()
	axarr[1].legend()
	plt.xlim(xmax=30)
	plt.show()

# shifts5
f, axarr = plt.subplots(2,sharex=True)
for j in [50]:
	allcondmax = []
	all_p_mat = []
	for g in [11.5,12.,12.5,13.]:#np.arange(7,12):
		for q in [11.5,12.,12.5,13.]:#np.arange(7,12):
			(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
			#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
			axarr[0].set_title('species m for mmax=50')
			axarr[1].plot(np.sum(p_mat_n_m[0],axis=1),label='g='+str(g)+',q='+str(q))
			axarr[1].set_title('species n for nmax=50')
	axarr[0].legend()
	axarr[1].legend()
	plt.xlim(xmax=30)
	plt.show()

# shifts6
f, axarr = plt.subplots(2,sharex=True)
for j in [50]:
	allcondmax = []
	all_p_mat = []
	for g in [9.5,10.,10.5,11.]:#np.arange(7,12):
		for q in [11.5,12.,12.5,13.]:#np.arange(7,12):
			(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
			#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
			axarr[0].set_title('species m for mmax=50')
			axarr[1].plot(np.sum(p_mat_n_m[0],axis=1),label='g='+str(g)+',q='+str(q))
			axarr[1].set_title('species n for nmax=50')
	axarr[0].legend()
	axarr[1].legend()
	plt.xlim(xmax=30)
	plt.show()
	
# shifts7
f, axarr = plt.subplots(2,sharex=True)
for j in [50]:
	allcondmax = []
	all_p_mat = []
	for g in [11.5,12.,12.5,13.]:#np.arange(7,12):
		for q in [9.5,10.,10.5,11.]:#np.arange(7,12):
			(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
			#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
			axarr[0].set_title('species m for mmax=50')
			axarr[1].plot(np.sum(p_mat_n_m[0],axis=1),label='g='+str(g)+',q='+str(q))
			axarr[1].set_title('species n for nmax=50')
	axarr[0].legend()
	axarr[1].legend()
	plt.xlim(xmax=30)
	plt.show()
	
# shifts8
f, axarr = plt.subplots(2,sharex=True)
for j in [50]:
	allcondmax = []
	all_p_mat = []
	for g in [11.5,12.,12.5,13.]:#np.arange(7,12):
		for q in [7.5,8.,8.5,9.]:#np.arange(7,12):
			(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
			#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
			axarr[0].set_title('species m for mmax=50')
			axarr[1].plot(np.sum(p_mat_n_m[0],axis=1),label='g='+str(g)+',q='+str(q))
			axarr[1].set_title('species n for nmax=50')
	axarr[0].legend()
	axarr[1].legend()
	plt.xlim(xmax=30)
	plt.show()

# shifts9
f, axarr = plt.subplots(2,sharex=True)
for j in [50]:
	allcondmax = []
	all_p_mat = []
	for g in [7.5,8.,8.5,9.]:#np.arange(7,12):
		for q in [11.5,12.,12.5,13.]:#np.arange(7,12):
			(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
			#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
			axarr[0].set_title('species m for mmax=50')
			axarr[1].plot(np.sum(p_mat_n_m[0],axis=1),label='g='+str(g)+',q='+str(q))
			axarr[1].set_title('species n for nmax=50')
	axarr[0].legend()
	axarr[1].legend()
	plt.xlim(xmax=30)
	plt.show()

#load Gillespie to compare:
p_mat = np.load('../../Gillespie/GeneRegulation_Gillespie/Data/p_n_m_1trajectory_endtime100000.npy')

f, axarr = plt.subplots(2,sharex=True)
for j in [50]:
	allcondmax = []
	all_p_mat = []
	for g in [7.254,11.5,12.5]:#np.arange(7,12):
		for q in [7.43,11.5,12.5]:#np.arange(7,12):
			(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
			#np.save("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy",p_mat_n_m[0])
			allcondmax.append(cond)
			all_p_mat.append(p_mat_n_m)
			axarr[0].plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
			axarr[0].set_title('species m for mmax=50 with shifts g='+str(g)+'and q='+str(q))
			axarr[1].plot(np.sum(p_mat_n_m[0],axis=1),label='g='+str(g)+',q='+str(q))
			axarr[1].set_title('species n for nmax=50 with shifts g='+str(g)+'and q='+str(q))
	axarr[0].legend()
	axarr[0].plot(dist_m[0],label='Gillespie')
	#axarr[0].errorbar(dist_m[0],dist_m[1],yerr=1.96/np.sqrt(num_traject)*std_m[1])
	axarr[1].legend()
	axarr[1].plot(dist_n[0],label='Gillespie')
	#axarr[1].errorbar(dist_n[0],dist_n[1],yerr=1.96/np.sqrt(num_traject)*std_n[1])
	plt.xlim(xmax=30)
	plt.show()

#def main_examine_shifts(nmax,mmax,jmax,kmax,g_mean,q_mean):
	## vector of initial marginal distribution p(n)
	#p_vec = p_vector(nmax)

	### vectors of creation rates g(n) and q(n) to find the mean rates for the shifts -> why needed?
	##q_vec = rate_vec(nmax,q_n0)
	##g_vec = rate_vec(nmax,g_n0)

	##q_mean = np.dot(q_vec,p_vec) # constant part of the creation rate of species m
	##g_mean = np.dot(g_vec,p_vec) # special in this example: creation rate of both species are regulated by species n with the same functional dependency
	#delta_q = delta_matrix(q_n0,nmax,q_mean) # diagonal matrix with entries (q_mean - q(n)) on the n^th diagonal element
	#delta_g = delta_matrix(g_n0,nmax,g_mean)
	#dual_base_mat_n = dual_mat(g_mean,nmax,jmax) # contains elements <j|n> TODO: not normalized! 
	#base_mat_n = mat(g_mean,nmax,jmax) # contains elements <n|j> - normalized :-)
	#base_mat_m = mat(q_mean,mmax,kmax) # contains elements <m|k> - normalized :-) 
	#new_delta_q = np.dot(np.dot(dual_base_mat_n.T,delta_q),base_mat_n) # delta_q in new base {|j>}
	#new_delta_g = np.dot(np.dot(dual_base_mat_n.T,delta_g),base_mat_n)
	#sub = subdiag(jmax,jmax)

	## gen will contain the coefficients of the generating function in the |j,m> base
	#gen = np.zeros((jmax,kmax))

	## first generating function: vector over indecies j with k=0
	#gen[:,0] = np.dot(dual_base_mat_n.T, p_vec)

	#for k in np.arange(1,kmax):
		#diag = np.zeros((jmax,jmax))
		#for j in np.arange(jmax):
			#diag[j,j] = j + rho*k
		#gen[:,k] = la.solve(-rho*(diag + np.dot(sub,new_delta_g)), np.dot(new_delta_q,gen[:,k-1]))

	#p_mat = np.dot(base_mat_n,np.dot(gen,base_mat_m.T)) # matrix containing the full distribution, i.e. for row n and column m it contains p(n,m)
	#return p_mat
	
##for j in [20,30,40,50,60]:
	##for g in np.arange(9,13):
		##for q in np.arange(9,13):
			##p_mat = main_examine_shifts(j,j,j,j,g,q)
			##np.save('Data/initial_hill2/p_'+str(j)+'_'+str(j)+'_'+str(j)+'_'+str(j)+'_shifts_'+str(g)+'_'+str(q)+'.npy',p_mat)
