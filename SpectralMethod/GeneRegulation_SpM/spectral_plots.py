import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pylab as p
from initializations import *
from functions import *
from main import *


allcond = []
all_p_mat = []
plt.figure()
for g in np.arange(6,13):
	for q in np.arange(6,13):
		(p_mat_n_m,cond) = ((spec_meth(g,q)[0],g,q),(max(spec_meth(g,q)[1]),g,q))
		allcond.append(cond)
		all_p_mat.append(p_mat_n_m)
		plt.title('marginal probability density of output proteins m')
		plt.xlabel('m')
		plt.ylabel('p(m)')
		plt.plot(np.sum(p_mat_n_m[0],axis=0),label='g='+str(g)+',q='+str(q))
plt.legend()
plt.show()

nmax = 30
mmax = 30
plt.figure()
plt.title('marginal probability density of input proteins n')
plt.xlabel('n')
plt.ylabel('p(n)')
for j in [20,30]:
	for g in np.arange(9,13):
		for q in np.arange(9,13):
			p_mat = np.load('Data/p_'+str(j)+'_'+str(j)+'_'+str(j)+'_'+str(j)+'_shifts_'+str(g)+'_'+str(q)+'.npy')
			p_n = np.sum(p_mat, axis=1) # sum of all elements in a row n gives value of the marginal distribution p(n)
			p_m = np.sum(p_mat, axis=0)
			plt.plot(np.sum(p_mat,axis=1),label='q='+str(q)+',g='+str(g)+'N,M,K,L='+str(j))
plt.legend()
plt.show()

plt.figure()
plt.title('marginal probability density of input proteins m')
plt.xlabel('m')
plt.ylabel('p(m)')
for j in [20,30]:
	for g in 2*np.arange(5,6):#[5,8,11,14]:
		for q in 2*np.arange(5,6):#[5,8,11,14]:
			p_mat = np.load("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(q)+".npy")
			p_n = np.sum(p_mat, axis=1) # sum of all elements in a row n gives value of the marginal distribution p(n)
			p_m = np.sum(p_mat, axis=0)
			plt.plot(np.sum(p_mat,axis=0),label='q='+str(q)+',g='+str(g))
plt.legend()
plt.show()

f, axarr = plt.subplots(2, 2)
for j in [50]:
	for g in np.arange(9,13):
		p_mat = np.load("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(7)+".npy")
		axarr[0,0].plot(np.sum(p_mat,axis=1),label='g='+str(g) )
		axarr[0,0].set_title('q=7')
		p_mat = np.load("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(9)+".npy")
		axarr[0,1].plot(np.sum(p_mat,axis=1),label='g='+str(g) )
		axarr[0,1].set_title('q=9')
		p_mat = np.load("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(11)+".npy")
		axarr[1,0].plot(np.sum(p_mat,axis=1),label='g='+str(g) )
		axarr[1,0].set_title('q=11')
		p_mat = np.load("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(12)+".npy")
		axarr[1,1].plot(np.sum(p_mat,axis=1),label='g='+str(g) )
		axarr[1,1].set_title('q=12')
	plt.legend()
	plt.ylabel('marginal probability density p(n)')
	plt.xlabel('protein number n')
	plt.show()
	
f, axarr = plt.subplots(2, 2)
for j in [50]:
	for g in np.arange(9,13):
		p_mat = np.load("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(7)+".npy")
		axarr[0,0].plot(np.sum(p_mat,axis=0),label='g='+str(g) )
		axarr[0,0].set_title('q=7')
		p_mat = np.load("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(9)+".npy")
		axarr[0,1].plot(np.sum(p_mat,axis=0),label='g='+str(g) )
		axarr[0,1].set_title('q=9')
		p_mat = np.load("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(11)+".npy")
		axarr[1,0].plot(np.sum(p_mat,axis=0),label='g='+str(g) )
		axarr[1,0].set_title('q=11')
		p_mat = np.load("Data/initial_hill2/p_"+str(j)+"_"+str(j)+"_"+str(j)+"_"+str(j)+"_shifts_"+str(g)+"_"+str(12)+".npy")
		axarr[1,1].plot(np.sum(p_mat,axis=0),label='g='+str(g) )
		axarr[1,1].set_title('q=12')
plt.legend()
plt.ylabel('marginal probability density p(m)')
plt.xlabel('protein number m')
plt.show()

#for g in np.arange(9,12):
	#for q in np.arange(9,12):
		#p_mat = np.load("../../SpectralMethod/GeneRegulation_SpM/Data/p_n_m_"+str(g)+"_"+str(q)+".npy")
		#p_n = np.sum(p_mat, axis=1) # sum of all elements in a row n gives value of the marginal distribution p(n)
		#p_m = np.sum(p_mat, axis=0) # sum of all elements in a column m gives value of the marginal distribution p(m)

		#n_mean = np.dot(np.arange(nmax),p_n) # compare with theoretically expected value: <n> = <g(n)>/1 where 1 is the constant degradation rate here
		#m_mean = np.dot(np.arange(mmax),p_m)

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
		
x = np.arange(0,nmax)
X, Y = p.meshgrid(x, x)
fig=plt.figure()
ax = Axes3D(fig)
ax.contour3D(X,Y,p_mat.T)
ax.set_xlabel('n')
ax.set_ylabel('m')
ax.set_zlabel('p(n,m)')
plt.show()

x = np.arange(0,nmax)
X, Y = p.meshgrid(x, x)
fig=plt.figure()
ax = Axes3D(fig)
ax.plot_wireframe(X,Y,p_mat)
ax.set_xlabel('m')
ax.set_ylabel('n')
ax.set_zlabel('p(n,m)')
plt.show()

plt.figure()
plt.title('marginal probability density of input proteins n')
plt.xlabel('n')
plt.ylabel('p(n)')
plt.plot(np.sum(p_mat,axis=1))
plt.legend()
plt.show()

plt.figure()
plt.title('marginal probability density of input proteins n')
plt.xlabel('m')
plt.ylabel('p(m)')
plt.plot(np.arange(nmax),np.sum(p_mat,axis=0))
plt.show()