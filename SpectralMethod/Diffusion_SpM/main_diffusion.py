import numpy as np
import math
from initializations_diffusion import *
from functions_diffusion import *
import numpy.linalg as la

def diffusion_sm_shift(N,a_bar,a_degger_bar):
	# matrices to shift above and below diagonal and solve initial equation
	diag_n_count = diag_count(N)
	delta = np.eye(N)*a_bar - np.dot(diag_n_count,np.tri(N,N,-1)-np.tri(N,N,-2))
	delta_degger = np.eye(N)*a_degger_bar - np.tri(N,N,1)+np.tri(N,N,0)

	# store generating function coefficients in this matrix:
	gen = np.zeros((N,N))

	# two initial vectors for k2=0 and k2=1 ??? this initial condition does not necessarily work!
	#for n1 in np.arange(N):
		#gen[n1,0] = 1./2.**(N-1)*np.float(math.factorial(N-1))/(math.factorial(N-1-n1)*math.factorial(n1))
	gen[:,0] = np.ones(N)
	gen[:,1] = np.dot(delta,gen[:,0])

	# calculate full matrix gen, i.e. for row n1 and column k2 it contains gen(n1,k2) at steady state
	for k2 in np.arange(1,N-1):
		gen[:,k2+1] = la.solve((k2+1)*delta_degger , -np.dot( np.dot(delta_degger,delta) + np.eye(N)*k2 ,gen[:,k2]) - np.dot(delta,gen[:,k2-1]))

	coeff_mat = mat(a_bar*a_degger_bar,N).T

	p_mat = np.dot(gen,coeff_mat)
	p_n1 = np.sum(p_mat, axis=1) # sum of each row n gives value of the marginal distribution p(n)
	p_n2 = np.sum(p_mat, axis=0) # sum of each column m gives value of the marginal distribution p(m)

	n1_mean = np.dot(np.arange(N),p_n1) # compare with theoretically expected value: 
	n2_mean = np.dot(np.arange(N),p_n2)
	print((n1_mean,n2_mean,np.sum(p_mat)))
	return p_mat
  
def save_data_shifts(N,a_bar_min,a_bar_max,a_degger_bar_min,a_degger_bar_max):
	for a_bar in np.arange(a_bar_min,a_bar_max):
		for a_degger_bar in np.arange(a_degger_bar_min,a_degger_bar_max):
			p_mat = diffusion_sm_shift(N,a_bar,a_degger_bar)
			np.save("/Data/p_"+str(N)+"_"+"shifts_"+str(a_bar)+"_"+str(a_degger_bar)+".npy",p_mat)  

p_mat = diffusion_sm_shift(50,6,6)

p_n1 = np.sum(p_mat, axis=1) # sum of each row n gives value of the marginal distribution p(n)
p_n2 = np.sum(p_mat, axis=0) # sum of each column m gives value of the marginal distribution p(m)

n1_mean = np.dot(np.arange(N),p_n1) # compare with theoretically expected value: 
n2_mean = np.dot(np.arange(N),p_n2)

##-----------------------recursion trial without shifting ladder operators

def diffusion_recurrence(N):
	diag = diag_count(N)
	A1 = np.tri(N,N,1)-np.tri(N,N,0)
	A2 = np.tri(N,N,2)-np.tri(N,N,1)

	p_mat_old = np.zeros((N,N))

	# two initial vectors:
	p_mat_old[:,0] = np.zeros(N)
	p_mat_old[N-1,0] = 1.
	p_mat_old[:,1] = np.dot(np.dot(A1,diag),p_mat_old[:,0])

	# calculate matrix containing the full distribution, i.e. for row n1 and column n2 it contains p(n1,n2) at steady state
	for n2 in np.arange(2,N):
		p_mat_old[:,n2] = 1./n2*(-np.dot(np.dot(A2,diag),p_mat_old[:,n2-2]) + (n2-1)*np.dot(A1,p_mat_old[:,n2-1]) + np.dot(np.dot(A1,diag),p_mat_old[:,n2-1]))

	p_mat_old = p_mat_old/np.sum(p_mat_old)
	return p_mat_old

#p_mat_old = diffusion_recurrence(N)

#p_n1_old = np.sum(p_mat_old, axis=1) # sum of each row n gives value of the marginal distribution p(n)
#p_n2_old = np.sum(p_mat_old, axis=0) # sum of each column m gives value of the marginal distribution p(m)

#n1_mean_old = np.dot(np.arange(N),p_n1_old) # compare with theoretically expected value: 
#n2_mean_old = np.dot(np.arange(N),p_n2_old)