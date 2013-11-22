import numpy as np
import math
from initializations import *
from functions import *
import numpy.linalg as la

def spec_meth(shift_g,shift_q):
	q_mean = shift_q # constant part of the creation rate of species m
	g_mean = shift_g # special in this example: creation rate of both species are regulated by species n with the same functional dependency
	delta_q = delta_matrix(q_n0,nmax,q_mean) # diagonal matrix with entries (q_mean - q(n)) on the n^th diagonal element
	delta_g = delta_matrix(g_n0,nmax,g_mean)
	dual_base_mat_n = dual_mat(g_mean,nmax,jmax) # contains elements <j|n> TODO: not normalized! 
	base_mat_n = mat(g_mean,nmax,jmax) # contains elements <n|j> - normalized :-)
	base_mat_m = mat(q_mean,mmax,kmax) # contains elements <m|k> - normalized :-) 
	new_delta_q = np.dot(np.dot(dual_base_mat_n.T,delta_q),base_mat_n) # delta_q in new base {|j>}
	new_delta_g = np.dot(np.dot(dual_base_mat_n.T,delta_g),base_mat_n)
	sub = subdiag(jmax,jmax)

	# gen will contain the coefficients of the generating function in the |j,m> base
	gen = np.zeros((jmax,kmax))

	# first generating function: vector over indecies j with k=0
	gen[:,0] = np.dot(dual_base_mat_n.T, p_vec)
	allcond = []
	for k in np.arange(1,kmax):
		diag = np.zeros((jmax,jmax))
		for j in np.arange(jmax):
			diag[j,j] = j + rho*k
		gen[:,k] = la.solve(-rho*(diag + np.dot(sub,new_delta_g)), np.dot(new_delta_q,gen[:,k-1]))
		eigval = la.eig(-rho*(diag + np.dot(sub,new_delta_g)))
		cond = np.float(np.max(np.abs(eigval[0])))/np.min(np.abs(eigval[0]))
		allcond.append(cond)
	p_mat = np.dot(base_mat_n,np.dot(gen,base_mat_m.T)) # matrix containing the full distribution, i.e. for row n and column m it contains p(n,m)
	return (p_mat,allcond)

#examine the condition of inverted matrices
allcond = []
for g in np.arange(6,13):
	for q in np.arange(6,13):
		cond = (max(spec_meth(g,q)[1]),g,q)
		allcond.append(cond)

def save_data_g_q():
	for g in np.arange(5,15):
		for q in np.arange(5,15):
			np.save("/Data/initial_hill2/p_"+str(nmax)+"_"+str(mmax)+"_"+str(jmax)+"_"+"kmax_"+"shifts_"+str(g)+"_"+str(q)+".npy",spec_meth(g,q))