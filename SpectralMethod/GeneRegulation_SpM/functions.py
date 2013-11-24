import numpy as np
import math
from initializations import *

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
	
def q_n0(n):
	""" auxiliary function for q with one argument"""
	return q_hill(n,n0)

def rate_vec(nmax,rate_func):
	""" mean value over the input protein number of the creation rate for species m, rate q(n) """
	rate_vec = np.zeros(nmax)
	for n in np.arange(nmax):
		rate_vec[n] = rate_func(n)
	return rate_vec

def g_n0(n):
	""" auxiliary function for g with one argument. For this example species n regulates itself with the same rate as species m is regulated by n"""
	return 7.# q_hill(n,n0)#

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
	for n in np.arange(1,nmax):
		p_vec[n] = p_initial(n)
	return p_vec/np.sum(p_vec)

#----------- matrices with coefficients between new {|j>} and old {|n>} base, recursion used in PNAS according to A13 and A14 in review

def mat(const,nmax,jmax):
	""" according to A13 find matrix with components ((<n=0|j=0>,<n=0|j=1>,<n=0|j=2>,...),(<n=1|j=0>,<n=1|j=1>,<n=1|j=2>,...),(<n=2|j=0>,<n=2|j=1>,<n=2|j=2>,...),...) n->row, j->column"""
	mat = np.zeros((nmax,jmax))
	for j in np.arange(jmax):
		mat[0,j] = np.exp(-const)*(-1)**j
		mat[1,j] = np.exp(-const)*(-1)**j*np.float(const-j)
		for n in np.arange(1,nmax-1):
				mat[n+1,j] = 1./np.float(n+1)*(np.float(const+n-j)*mat[n,j] - np.float(const)*mat[n-1,j])
	return mat

def dual_mat(const,nmax,jmax):
	""" according to A14 find matrix with components ((<j=0|n=0>,<j=0|n=1>,<j=0|n=2>,...),(<j=1|n=0>,<j=1|n=1>,<j=1|n=2>,...),...) j->row, n->column"""
	mat = np.zeros((jmax,nmax))
	for j in np.arange(jmax):
		mat[j,0] = (-const)**j/math.factorial(j)
		mat[j,1] = (-const)**j*((np.float(const)-np.float(j))/np.float(const))/math.factorial(j)# <j=0|n> = 1 for all n
		for n in np.arange(1,nmax-1):
			mat[j,n+1] = 1./np.float(const)*(np.float(const+n-j)*mat[j,n]-np.float(n)*mat[j,n-1])
	return mat

#????? according to Mugler, Walczak review, this is less numerically stable:

#def mat(constant,nmax,jmax):
	#""" according to A10 find matrix with components ((<n=0|j=0>,<n=0|j=1>,<n=0|j=2>,...),(<n=1|j=0>,<n=1|j=1>,<n=1|j=2>,...),(<n=2|j=0>,<n=2|j=1>,<n=2|j=2>,...),...) n->row, j->column"""
	#mat = np.zeros((nmax,jmax))
	#for j in np.arange(jmax):
		#mat[0,j] = np.exp(-constant)*(-1)**j
		#for n in np.arange(nmax-1):
			#mat[n+1,j] = np.float(j)/np.float(n+1)*mat[n,j-1] + np.float(constant)/np.float(n+1)*mat[n,j] 
	#return mat

#def dual_mat(const,nmax,jmax):
	#""" according to A11 find matrix with components ((<j=0|n=0>,<j=0|n=1>,<j=0|n=2>,...),(<j=1|n=0>,<j=1|n=1>,<j=1|n=2>,...),...) j->row, n->column"""
	#mat = np.zeros((jmax,nmax))
	#for j in np.arange(jmax):
		#mat[j,0] = (-const)**j/math.factorial(j)
		#for n in np.arange(0,nmax-1):
			#if j==0:
				#mat[j,n+1] = mat[j,n]
			#else:
				#mat[j,n+1] = mat[j-1,n] + mat[j,n]
	#return mat

def delta_matrix(q_func, nmax, q_bar):
	""" returns the deviation matrix in the {|n>} base which has only diagonal entries (q_bar-q_n0(n)) in the |n> base """
	delta = np.identity(nmax)
	for i in np.arange(nmax):
		delta[i][i] = q_func(i) - q_bar
	return delta

def subdiag(nmax,jmax):
	""" returns a matrix of shape (nmax,jmax) with ones on the subdiagonal and zeros elsewhere"""
	return np.tri(nmax,jmax,-1)-np.tri(nmax,jmax,-2)



#----------- algebra with one type of recursion - didn't converge!?

def input_base_coeff_n0(j): 
	""" calculate <0|j> -> will be a different value for another (correct hilbert space?) scalar product! TODO """
	return (-1)**j*np.exp(-g)
	
def recursion_rel(coeff_minus_one, coeff_minus_two, n, j):
	""" recursion relation for the coefficients <n|j> """
	if(n==0):
		return input_base_coeff_n0(j)
	else:
		return (g - j + n - 1)/n*coeff_minus_one - g/n*coeff_minus_two
	
def input_base_coeff_vec(nmax,j):
	""" calculate a vector of length nmax containing <n|j> in each row for fix j """
	vec = np.zeros(nmax)
	vec[0] = input_base_coeff_n0(j)
	vec[1] = vec[0]*(g-j)
	for n in np.arange(2,nmax):
		vec[n] = (g - j + n - 1)/n*vec[n-1] - g/n*vec[n-2]
	return vec

def input_dual_coeff_n0(j): 
	""" calculate <j|0> -> will be a different value for another (correct hilbert space?) scalar product! TODO """
	if(j==0):
		return 1.
	else:
		return (-g)**j/math.factorial(j)

def input_dual_vec(nmax,j):
	""" calculate a vector of length nmax containing <j|n> in each row for fix j """
	vec = np.zeros(nmax)
	vec[0] = input_dual_coeff_n0(j)
	vec[1] = vec[0]*(1-j/g)
	for n in np.arange(2,nmax):
		vec[n] = (g - j + n - 1)/g*vec[n-1] - n/g*vec[n-2]
	return vec

def output_base_coeff_m0(k):
	""" calculate <m=0|k> -> might be a different value for another (correct hilbert space?) scalar product??? TODO """
	q = q_bar(nmax) 
	return (-1)**k*np.exp(-q)
	
def output_base_coeff_vec(mmax,k):
	""" calculate a vector of length mmax containing <m|k> in each row for fix k """
	q = q_bar(nmax)
	vec = np.zeros(mmax)
	vec[0] = output_base_coeff_m0(k)
	vec[1] = vec[0]*(q-k)
	for m in np.arange(2,mmax):
		vec[m] = (q - k + m - 1)/m*vec[m-1] - q/m*vec[m-2]
	return vec

def vec_matrix(jmax, nmax, vec_func):
	""" returns a matrix of shape (nmax,jmax). Each column is a vector of length nmax, where each vector is created with the function vec_func(nmax,j) """
	v_mat = np.zeros((nmax,jmax))
	for j in np.arange(jmax):
		v_mat[:][j] = vec_func(nmax,j)
	return v_mat
