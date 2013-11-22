import numpy as np
from initializations_cascade import *
from functions import *
from main_GeneReg_cascade import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pylab as p

#n1max = 40
#n2max = 40
#n3max = 40
p_n3_n1_n2 = p_mat_tensor #p_n3_n1_n2 = np.load('p_mat_'+str(n1max)+'_'+str(n2max)+'_'+str(n3max)+'_35_35_35.npy')
p_n3_n1 = np.sum(p_n3_n1_n2,axis=2) # marginal distribution p(n3,n1)
p_n3_n2 = np.sum(p_n3_n1_n2,axis=1) # marginal distribution p(n3,n2)
p_n1_n2 = np.sum(p_n3_n1_n2,axis=0)
p_n1 = np.sum(p_n3_n1,axis=0)
p_n3 = np.sum(p_n3_n1,axis=1)
p_n2 = np.sum(p_n3_n2,axis=0)

x = np.arange(0,n1max)
y = np.arange(0,n2max)
X, Y = p.meshgrid(x, x)
fig=plt.figure()
ax = Axes3D(fig)
ax.plot_wireframe(X,Y,p_n1_n2[:,:n1max])
ax.set_xlabel('n1')
ax.set_ylabel('n2')
ax.set_zlabel('p(n1,n2)')
p.show()

x = np.arange(0,n1max)
y = np.arange(0,n1max)
X, Y = p.meshgrid(x, y)
fig=plt.figure()
ax = Axes3D(fig)
ax.plot_wireframe(X,Y,p_n3_n1[:n1max,:])
ax.set_xlabel('n3')
ax.set_ylabel('n1')
ax.set_zlabel('p(n3,n1)')
p.show()

plt.figure()
plt.title('marginal probability density of input proteins n1')
plt.xlabel('n1')
plt.ylabel('p(n1)')
plt.plot(np.arange(n1max),p_n1)
plt.plot(np.arange(n1max),p_vec,'r')
plt.show()

plt.figure()
plt.title('marginal probability density of output protein n3')
plt.xlabel('n3')
plt.ylabel('p(n3)')
plt.plot(np.arange(n3max),p_n3)
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