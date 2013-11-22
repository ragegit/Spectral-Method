import numpy as np
from initializations_diffusion import *
from functions_diffusion import *
from main_diffusion import *
from main_diffusion_change_kmax import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pylab as p

#x = np.arange(0,N)
#X, Y = p.meshgrid(x, x)
#fig=plt.figure()
#ax = Axes3D(fig)
#ax.contour3D(X,Y,p_mat.T)
#ax.set_xlabel('n1')
#ax.set_ylabel('n2')
#ax.set_zlabel('p(n1,n2)')
#p.show()

def plotting(N,K,a_bar,a_degger_bar):
	p_mat = diffusion_sm_shift_K(N,K,a_bar,a_degger_bar)
	x = np.arange(0,N)
	#X, Y = p.meshgrid(x, x)
	#fig=plt.figure()
	#plt.ylim(-0.8,1.5)
	#plt.title('full distribution from recurrence equation after base switch')
	#ax = Axes3D(fig)
	#ax.plot_wireframe(X,Y,p_mat)
	#ax.set_xlabel('n1')
	#ax.set_ylabel('n2')
	#ax.set_zlabel('p(n1,n2)')
	#p.show()

	plt.figure()
	plt.ylim(-0.01,0.5)
	plt.title('marginal probability density at site 1')
	plt.xlabel('n1')
	plt.ylabel('p(n1)')
	x = np.arange(N)
	plt.plot(x,map(Bernoulli_dist,x,np.ones(N)*(N-1)),'b')
	plt.plot(np.sum(p_mat,axis=0),'r')
	plt.show()

	plt.figure()
	plt.ylim(-0.01,0.5)
	plt.title('marginal probability density at site 2')
	plt.xlabel('n2')
	plt.ylabel('p(n2)')
	plt.plot(x,map(Bernoulli_dist,x,np.ones(N)*(N-1)),'b')
	plt.plot(np.sum(p_mat,axis=1),'r')
	plt.show()
	
def plot_data_shifts(N,K,a_bar_min,a_bar_max,a_degger_bar_min,a_degger_bar_max):
	for a_bar in np.arange(a_bar_min,a_bar_max):
		for a_degger_bar in np.arange(a_degger_bar_min,a_degger_bar_max):
			p_mat = np.load("Data/p_"+str(N)+"_"+str(K)+"_"+"shifts_"+str(a_bar)+"_"+str(a_degger_bar)+".npy")
			plt.figure()
			plt.ylim(-0.01,0.5)
			plt.title('marginal probability density at site 1')
			plt.xlabel('n1')
			plt.ylabel('p(n1)')
			x = np.arange(N)
			plt.plot(x,map(Bernoulli_dist,x,np.ones(N)*(N-1)),'b')
			plt.plot(np.arange(N),np.sum(p_mat,axis=0),'r')
			plt.show()

			plt.figure()
			plt.ylim(-0.01,0.5)
			plt.title('marginal probability density at site 2')
			plt.xlabel('n2')
			plt.ylabel('p(n2)')
			plt.plot(x,map(Bernoulli_dist,x,np.ones(N)*(N-1)),'b')
			plt.plot(np.arange(N),np.sum(p_mat,axis=1),'r')
			plt.show()
	
plotting(80,100,5.5,5.5)

### ---------------- plot the direct recurrence results
x = np.arange(0,N)
X, Y = p.meshgrid(x, x)
fig=plt.figure()
plt.title('full distribution from direct recurrence equation')
ax = Axes3D(fig)
ax.plot_wireframe(X,Y,p_mat_old)
ax.set_xlabel('n1')
ax.set_ylabel('n2')
ax.set_zlabel('p(n1,n2)')
p.show()

plt.figure()
plt.title('marginal probability density at site 1 with direct recurrence')
plt.xlabel('n1')
plt.ylabel('p(n1)')
x = np.arange(N)
plt.plot(x,np.sum(p_mat_old,axis=0),label='Spectral Method')
#plt.plot(x,np.abs(np.sum(p_mat_old,axis=0))/np.sum(np.abs(np.sum(p_mat_old,axis=0))))
plt.plot(x,map(Bernoulli_dist,x,np.ones(N)*(N-1)),label='analytical result')
# plot the analytically calculated distribution against the simulation result
plt.legend()
plt.show()

plt.figure()
plt.title('marginal probability density at site 2 with direct recurrence')
plt.xlabel('n2')
plt.ylabel('p(n2)')
plt.plot(np.arange(N),np.sum(p_mat_old,axis=1),label='Spectral Method')
plt.legend()
plt.show()