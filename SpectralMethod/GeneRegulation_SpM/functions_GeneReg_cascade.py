import numpy as np
import math
from initializations_cascade import *

def hill(n,nu,n0):
	""" hill function """
	return np.float(n**nu)/(n**nu + n0**nu)

# creation rates------------------------------------

def q_step(n,n0):
	""" n dependend creation rate for output species """
	#return qminus + (qplus - qminus)*hill(n,nu,n0)
	if n<=n0:
		return qminus
	else:
		return qplus
	
def q_hill(n,n0):
	""" n dependend creation rate for output species """
	#return qminus + (qplus - qminus)*hill(n,nu,n0)
	return np.float(qminus*n0**nu + qplus*n**nu)/(n**nu + n0**nu)
	
def q2(n):
	""" creation rate of species 2 depending on the number of particles of species 1"""
	return q_hill(n,n0)

def q3(n):
	""" creation rate of species 3 depending on the number of particles of species 2"""
	return q_hill(n,n0)


def rate_vec(nmax,rate_func):
	""" mean value over the input protein number of the creation rate for species 2 or 3 with rate q2(n1) or q3(n2)"""
	rate_vec = np.zeros(nmax)
	for n in np.arange(nmax):
		rate_vec[n] = rate_func(n)
	return rate_vec

def g_n0(n):
	""" auxiliary function for g with one argument. For this example species n regulates itself with the same rate as species m is regulated by n"""
	return 8.

def product_g(n):
	""" auxiliary function to create a factorial for a function via it's natural number argument: g(0)*g(1)*g(2)*...*g(n)"""
	if(n>0):
		return product_g(n-1)*g_n0(n)
	else:
		return g_n0(0)
	
def facto(n):
	""" auxiliary function: calculate factorial recursively"""
	if(n>0):
		return facto(n-1)*n
	else:
		return 1
# probabilities---------------------------

def p_initial(n):
	""" steady state distribution for one independent species with autoregulation with rate g(n) - normalized in full distribution vecor below"""
	return 1./math.factorial(n)*product_g(n-1)

def p_vector(nmax):
	""" vector with input probabilities p(n) in each entry """
	p_vec = np.zeros(nmax)
	for n in np.arange(nmax):
		p_vec[n] = p_initial(n)
	return p_vec/np.sum(p_vec)

#----------- matrices with coefficients between new {|j>} and old {|n>} base, recursion used in PNAS

def mat(const,nmax,jmax):
	""" find matrix with components ((<n=0|j=0>,<n=0|j=1>,<n=0|j=2>,...),(<n=1|j=0>,<n=1|j=1>,<n=1|j=2>,...),(<n=2|j=0>,<n=2|j=1>,<n=2|j=2>,...),...) n->row, j->column"""
	mat = np.zeros((nmax,jmax))
	for n in np.arange(nmax):
		mat[n,0] = np.exp(-const)*const**n/math.factorial(n)
		for j in np.arange(jmax-1):
			if(n>0):
				mat[n,j+1]=mat[n-1,j]-mat[n-1,j]
			else:
				mat[n,j+1]=(-1)**(j+1)*np.exp(-const)
	return mat

def dual_mat(const,nmax,jmax):
	""" find matrix with components ((<j=0|n=0>,<j=1|n=0>,<j=2|n=0>,...),(<j=0|n=1>,<j=1|n=1>,<j=2|n=1>,...),...) n->row, j->column"""
	mat = np.zeros((nmax,jmax))
	for n in np.arange(nmax):
		mat[n,0] = 1 # <j=0|n> = 1 for all n
		for j in np.arange(jmax-1):
			if(n>0):
				mat[n,j+1] = mat[n-1,j] + mat[n-1,j+1]
			else:
				mat[n,j+1] = (-const)**(j+1)/math.factorial(j+1)
	return mat

def delta_matrix(q_func, nmax, q_bar):
	""" returns the deviation matrix in the {|n>} base which has only diagonal entries (q_bar-q_n0(n)) in the |n> base """
	delta = np.identity(nmax)
	for i in np.arange(nmax):
		delta[i][i] = q_func(i) - q_bar
	return delta

def subdiag(nmax,jmax):
	""" returns a matrix of shape (nmax,jmax) with ones on the subdiagonal and zeros elsewhere"""
	return np.tri(nmax,jmax,-1)-np.tri(nmax,jmax,-2)

