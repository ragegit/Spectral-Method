import numpy as np
import math
from initializations_diffusion import *

def subdiag(nmax,jmax):
	""" returns a matrix of shape (nmax,jmax) with ones above the diagonal and zeros elsewhere"""
	return np.tri(nmax,jmax,1)-np.tri(nmax,jmax,0)

def diag_count(n):
	""" returns a matrix of shape (n,n) with n at diagonal entry n and zeros elsewhere """
	diag = np.zeros((n,n))
	for j in np.arange(n):
		diag[j,j] = j
	return diag

def diag_k(a_bar, a_degger_bar,k,n):
	""" returns a matrix of shape (n,n) with (a_bar*a_degger+k) at diagonal entry n and zeros elsewhere"""
	diag = (a_bar*a_degger_bar+k)*np.eye(n)
	return diag

#----------- matrices with coefficients between new {|k2>} and old {|n2>} base, recursion used in PNAS

def mat(const,N):
	""" find matrix with components ((<n2=0|k2=0>,<n2=0|k2=1>,<n2=0|k2=2>,...),(<n2=1|k2=0>,<n2=1|k2=1>,<n2=1|k2=2>,...),(<n2=2|k2=0>,<n2=2|k2=1>,<n2=2|k2=2>,...),...) n->row, j->column"""
	mat = np.zeros((N,N))
	for n in np.arange(N):
		mat[n,0] = np.exp(-const)*const**n/math.factorial(n)
		for k in np.arange(N-1):
			if(n>0):
				mat[n,k+1]=mat[n-1,k]-mat[n-1,k]
			else:
				mat[n,k+1]=(-1)**(k+1)*np.exp(-const)
	return mat

def mat_K(const,N,K):
	""" find matrix with components ((<n2=0|k2=0>,<n2=0|k2=1>,<n2=0|k2=2>,...),(<n2=1|k2=0>,<n2=1|k2=1>,<n2=1|k2=2>,...),(<n2=2|k2=0>,<n2=2|k2=1>,<n2=2|k2=2>,...),...) n->row, j->column"""
	mat = np.zeros((N,K))
	for n in np.arange(N):
		mat[n,0] = np.exp(-const)*const**n/math.factorial(n)
		for k in np.arange(K-1):
			if(n>0):
				mat[n,k+1]=mat[n-1,k]-mat[n-1,k]
			else:
				mat[n,k+1]=(-1)**(k+1)*np.exp(-const)
	return mat

# compare to analytical results
def Bernoulli_dist(k,N):
	"""function to dray the distribution"""
	return 1./2.**N*math.factorial(N)/math.factorial(N-k)*1./math.factorial(k)
