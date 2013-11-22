import numpy as np
import math
from initializations_cascade import *
from functions_GeneReg_cascade import *
from initializations import *
from functions import *
import numpy.linalg as la
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pylab as p

#------------------------------------------------two species:
# at first calculate the marginal distribution from a two species cascade in order to initialize the generating function of the three species cascade:

# vector of initial marginal distribution p(n) - check if this marginal dist. is the same as p_n1
p_vec = p_vector(n1max)
# vector of constant creation rate of species 1 and non constant creation rates of species 2
g_vec = rate_vec(n1max,g_n0)
q2_vec = rate_vec(n1max,q2)
# constant part of the creation rate of species n2
q2_mean = np.dot(q2_vec,p_vec) # constant part of the creation rate of species n2
g_mean = np.dot(g_vec,p_vec)
# diagonal matrix with entries (q2_mean - q2(n1)) on the n1^th diagonal element
delta_g = delta_matrix(g_n0,n1max,g_mean)
delta_q2_n1 = delta_matrix(q2,n1max,q2_mean)

dual_base_mat_n1 = dual_mat(g_mean,n1max,j1max) # contains elements <j1|n1> TODO: not normalized? 
base_mat_n1 = mat(g_mean,n1max,j1max) # contains elements <n1|j1> - normalized :-)
base_mat_n2 = mat(q2_mean,n2max,j2max) # contains elements <n2|j2> - normalized :-) 
new_delta_q2_n1 = np.dot(np.dot(dual_base_mat_n1.T,delta_q2_n1),base_mat_n1)# delta_q in new base {|j>}
new_delta_g = np.dot(np.dot(dual_base_mat_n1.T,delta_g),base_mat_n1)
sub = subdiag(j1max,j1max)

# gen will contain the coefficients of the generating function in the |j1,j2> base
gen = np.zeros((j1max,j2max))

# first generating function: vector over indecies j1 with j2=0
gen[:,0] = np.dot(dual_base_mat_n1.T, p_vec)

diag = np.zeros((j1max,j1max))
for j2 in np.arange(1,j2max):
	for j1 in np.arange(j1max):
		diag[j1,j1] = j1 + rho2*j2
	gen[:,j2] = la.solve(-rho2*(diag + np.dot(sub,new_delta_g)), np.dot(new_delta_q2_n1,gen[:,j2-1]))
	
p_mat_n1_n2 = np.dot(base_mat_n1,np.dot(gen,base_mat_n2.T))# matrix containing the full distribution of a two species cascade, i.e. for row n and column m it contains p(n,m)
p_n1 = np.sum(p_mat_n1_n2, axis=1) # sum of each row n gives value of the marginal distribution p(n)
p_n2 = np.sum(p_mat_n1_n2, axis=0) # sum of each column m gives value of the marginal distribution p(m)

#x = np.arange(0,n1max)
#y = np.arange(0,n2max)
#X, Y = p.meshgrid(x, y)
#fig=plt.figure()
#ax = Axes3D(fig)
#ax.plot_wireframe(X,Y,p_mat_n1_n2)
#ax.set_xlabel('n2')
#ax.set_ylabel('n1')
#ax.set_zlabel('p(n1,n2)')
#p.show()

plt.figure()
plt.title('marginal probability density of input proteins n')
plt.xlabel('n')
plt.ylabel('p(n)')
plt.plot(np.arange(n2max),np.sum(p_mat_n1_n2,axis=0))
plt.plot(np.arange(n1max),p_vec)
plt.show()
# --------------------------------------- now use gen to initialize three species

# vector of non constant creation rates of species 3
q3_vec = rate_vec(n2max,q3)
q3_mean = np.dot(q3_vec,p_n2)
delta_q3_n2 = delta_matrix(q3,n2max,q3_mean)
dual_base_mat_n2 = dual_mat(q2_mean,n2max,j2max) # contains elements <j2|n2> TODO: not normalized! 
new_delta_q3_n2 = np.dot(np.dot(dual_base_mat_n2.T,delta_q3_n2),base_mat_n2)
base_mat_n3 = mat(q3_mean,n3max,j3max) # contains elements <n3|j3> - normalized :-) 

# matrices needed to solve recursion relation:
diag_j1 = np.zeros((j1max*j2max,j1max*j2max))
diag_j2 = np.zeros((j1max*j2max,j1max*j2max))
Delta_q2_j1 = np.zeros((j1max*j2max,j1max*j2max))
Delta_q3_j2 = np.zeros((j1max*j2max,j1max*j2max))
B_j2 = np.tri(j1max*j2max,j1max*j2max,-j1max-1)-np.tri(j1max*j2max,j1max*j2max,-j1max)
#Base_mat_n1 = np.zeros((n1max*j2max*j3max,j1max*j2max*j3max))
base1_help = np.zeros((n1max*j2max,j1max*j2max))
#Base_mat_n2 = np.zeros((n1max*n2max*j3max,n1max*j2max*j3max))
base2_help = np.zeros((n1max*n2max,n1max*j2max))
#Base_mat_n3 = base_mat_n3.T

for j in np.arange(j1max*j2max):
	diag_j1[j,j] = j%j1max
	diag_j2[j,j] = int(j/j1max)
		
for j in np.arange(j2max):
	Delta_q2_j1[j*j1max:(j+1)*j1max,j*j1max:(j+1)*j1max] = new_delta_q2_n1
	Delta_q3_j2[j*j1max:(j+1)*j1max,j*j1max:(j+1)*j1max] = new_delta_q3_n2[int(j/j2max),int(j/j2max)]*np.eye(j1max)
	base1_help[j*n1max:(j+1)*n1max,j*j1max:(j+1)*j1max] = base_mat_n1
	for i in np.arange(n2max):
		base2_help[i*n1max:(i+1)*n1max,j*n1max:(j+1)*n1max] = base_mat_n2[i%n2max,j%n2max]*np.eye(n1max)

#for j in np.arange(j3max):
	#Base_mat_n1[j*(n1max*j2max):(j+1)*(n1max*j2max),j*(j1max*j2max):(j+1)*(j1max*j2max)] = base1_help
	#Base_mat_n2[j*(n1max*n2max):(j+1)*(n1max*n2max),j*(n1max*j2max):(j+1)*(n1max*j2max)] = base2_help

# gen_cascade will contain the coefficients of the generating function in the |j1,j2,j3> base
gen_cascade = np.zeros((j1max*j2max,j3max))

# first generating function: vector over indecies j with k=0
gen_cascade[:,0] = gen.flatten('F')

for j3 in np.arange(1,j3max):
	gen_cascade[:,j3] = la.solve(-(diag_j1+rho2*diag_j2+rho3*j3*np.eye(j1max*j2max)+rho2*np.dot(Delta_q2_j1,B_j2)),rho3*np.dot(Delta_q3_j2,gen_cascade[:,j3-1]))

p_mat = np.dot(base_mat_n3,np.dot(base2_help,np.dot(base1_help,gen_cascade)).T)
p_mat_tensor = np.reshape(p_mat,(n3max,n1max,n2max))

p_n3_n1 = np.sum(p_mat_tensor,axis=2) # marginal distribution p(n3,n1)
p_n3_n2 = np.sum(p_mat_tensor,axis=1) # marginal distribution p(n3,n2)
p_n1_n2 = np.sum(p_mat_tensor,axis=0) # marginal distribution p(n3,n2)
p_n1 = np.sum(p_n3_n1,axis=0)
p_n3 = np.sum(p_n3_n1,axis=1)
p_n2 = np.sum(p_n3_n2,axis=0)

n1_mean = np.dot(np.arange(n1max),p_n1) # compare with theoretically expected value: <n> = <g(n)>/1 where 1 is the constant degradation rate here
n2_mean = np.dot(np.arange(n2max),p_n2)
n3_mean = np.dot(np.arange(n3max),p_n3)