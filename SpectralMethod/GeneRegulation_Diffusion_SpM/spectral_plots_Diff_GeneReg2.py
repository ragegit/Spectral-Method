import numpy as np
from initializations_Diff_GeneReg import *
from functions_Diff_GeneReg import *
from main_Diff_GeneReg import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pylab as p

p_mat_m2_n2_m1_n1 = p_mat_m2_n1_m1_n2/np.sum(p_mat_m2_n2_m1_n1)

p_n2_m1_n1 = np.sum(p_mat_m2_n2_m1_n1,axis = 0)
p_m1_n1 = np.sum(p_n2_m1_n1,axis = 0)
p_n1 = np.sum(p_m1_n1, axis = 0)
p_n2 = np.sum(np.sum(p_n2_m1_n1,axis = 1),axis=1)
p_m1 = np.sum(p_m1_n1,axis = 1)

x = np.arange(0,M)
y = np.arange(0,M)
X, Y = p.meshgrid(x, y)
fig=plt.figure()
ax = Axes3D(fig)
ax.plot_wireframe(X,Y,p_m1_n1[:,:])
ax.set_xlabel('m1')
ax.set_ylabel('n1')
ax.set_zlabel('p(m1,n1)')
p.show()

plt.figure()
plt.title('marginal probability density of input proteins n1')
plt.xlabel('n1')
plt.ylabel('p(n1)')
plt.plot(np.arange(N),p_n1)
plt.plot(np.arange(N),p_vec,'r')
plt.show()

plt.figure()
plt.title('marginal probability density of output protein n2')
plt.xlabel('n2')
plt.ylabel('p(n2)')
plt.plot(np.arange(N),p_n2)
plt.show()

plt.figure()
plt.title('marginal probability density of output protein m1')
plt.xlabel('m1')
plt.ylabel('p(m1)')
plt.plot(np.arange(M),p_m1)
plt.show()

#x = np.arange(0,nmax)
#X, Y = p.meshgrid(x, x)
#fig=plt.figure()
#ax = Axes3D(fig)
#ax.plot_wireframe(X,Y,p_mat)
#ax.set_xlabel('m')
#ax.set_ylabel('n')
#ax.set_zlabel('p(n,m)')
#p.show()

#plt.figure()
#plt.title('marginal probability density of input proteins n')
#plt.xlabel('n')
#plt.ylabel('p(n)')
#plt.plot(np.arange(nmax),np.sum(p_mat,axis=0))
#plt.show()

#plt.figure()
#plt.title('marginal probability density of input proteins n')
#plt.xlabel('m')
#plt.ylabel('p(m)')
#plt.plot(np.arange(nmax),np.sum(p_mat,axis=1))
#plt.show()