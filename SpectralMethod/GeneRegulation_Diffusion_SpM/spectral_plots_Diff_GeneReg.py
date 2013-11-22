#from initializations_Diff_GeneReg import *
#from functions_Diff_GeneReg import *
#from main_Diff_GeneReg import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import pylab as p

p_m2_n2_m1_n1 = np.load('Data/initial1/p_'+str(N)+'_'+str(M)+'_'+str(K)+'_'+str(L)+'_shifts_'+str(g)+'_'+str(q)+'.npy')
#p_m2_n2_m1_n1 = np.reshape(np.load('data_initial2/p_18_18_18_18_shifts_'+str(q)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2

p_n2_m1_n1 = np.sum(p_m2_n2_m1_n1,axis = 0)
p_m1_n1 = np.sum(p_n2_m1_n1,axis = 0)
p_n1 = np.sum(p_m1_n1, axis = 0)
p_n2 = np.sum(np.sum(p_n2_m1_n1,axis = 1),axis=1)
p_m1 = np.sum(p_m1_n1,axis = 1)

#p_n1_m1_n2 = np.sum(p_mat_m2_n1_m1_n2,axis = 0)
#p_m1_n2 = np.sum(p_n1_m1_n2,axis = 0)
##p_n1 = np.sum(np.sum(p_n1_m1_n2,axis = 1), axis = 1)
#p_n2 = np.sum(p_m1_n2,axis = 0)
#p_m1 = np.sum(p_m1_n2,axis = 1)
#p_n1 = np.sum(np.sum(np.sum(p_n1_n2_m1_m2_Spec, axis = 1), axis = 1), axis = 1)

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
for g in np.arange(9,12):
	#for q in np.arange(9,12):
	p_m2_n2_m1_n1 = np.save('Data/initial1/p_'+str(N)+'_'+str(M)+'_'+str(K)+'_'+str(L)+'_shifts_'+str(g)+'_'+str(q)+'.npy')
	p_m1_n1 = np.sum(p_n2_m1_n1,axis = 0)
	p_n1 = np.sum(p_m1_n1, axis = 0)
	plt.plot(p_n1,label='q='+str(q)+',g='+str(g) )
plt.legend()
plt.show()

f, axarr = plt.subplots(3, 1)
for g in np.arange(9,12):
	#for q in np.arange(9,12):
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(9)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[0].plot(p_n1,label='g='+str(g) )
	axarr[0].set_title('q=9')
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(10)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[1].plot(p_n1,label='g='+str(g) )
	axarr[1].set_title('q=10')
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(11)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[2].plot(p_n1,label='g='+str(g) )
	axarr[2].set_title('q=11')
plt.legend()
plt.ylabel('marginal probability density p(n1)')
plt.xlabel('protein number n1')
plt.show()

f, axarr = plt.subplots(3, 1)
for g in np.arange(9,12):
	#for q in np.arange(9,12):
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(9)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[0].plot(p_n2,label='g='+str(g) )
	axarr[0].set_title('q=9')
	plt.title('marginal distribution p(n2) versus n2')
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(10)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[1].plot(p_n2,label='g='+str(g) )
	axarr[1].set_title('q=10')
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(11)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[2].plot(p_n2,label='g='+str(g) )
	axarr[2].set_title('q=11')
plt.legend()
plt.ylabel('marginal probability density p(n2)')
plt.xlabel('protein number n2')
plt.show()

f, axarr = plt.subplots(3, 1)
for g in np.arange(9,12):
	#for q in np.arange(9,12):
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(9)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[0].plot(p_m1,label='g='+str(g) )
	axarr[0].set_title('q=9')
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(10)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[1].plot(p_m1,label='g='+str(g) )
	axarr[1].set_title('q=10')
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(11)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[2].plot(p_m1,label='g='+str(g) )
	axarr[2].set_title('q=11')
plt.legend()
plt.ylabel('marginal probability density p(m1)')
plt.xlabel('protein number m1')
plt.show()

f, axarr = plt.subplots(3, 1)
for g in np.arange(9,12):
	#for q in np.arange(9,12):
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(9)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[0].plot(p_m2,label='g='+str(g) )
	axarr[0].set_title('q=9')
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(10)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[1].plot(p_m2,label='g='+str(g) )
	axarr[1].set_title('q=10')
	p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(11)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
	p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
	p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
	p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
	p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
	p_m1 = np.sum(p_m2_m1, axis = 1)
	p_m2 = np.sum(p_m2_m1, axis = 0)
	axarr[2].plot(p_m2,label='g='+str(g) )
	axarr[2].set_title('q=11')
plt.legend()
plt.ylabel('marginal probability density p(m2)')
plt.xlabel('protein number m2')
plt.show()

plt.figure()
plt.title('marginal probability density of output protein m1')
plt.xlabel('m1')
plt.ylabel('p(m1)')
for g in np.arange(9,12):
	#for q in np.arange(9,12):
		p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(q)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
		p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
		p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
		p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
		p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
		p_m1 = np.sum(p_m2_m1, axis = 1)
		p_m2 = np.sum(p_m2_m1, axis = 0)
		plt.plot(p_m1,label='q='+str(q)+',g='+str(g) )
plt.legend()
plt.show()

plt.figure()
plt.title('marginal probability density of output protein m2')
plt.xlabel('m2')
plt.ylabel('p(m2)')
#for g in np.arange(9,12):
	for q in np.arange(9,12):
		p_m2_n1_m1_n2 = np.reshape(np.load('GeneRegulation_Diffusion_SpM/data_initial2/p_18_18_18_18_shifts_'+str(q)+'_'+str(g)+'.npy'),(18,18,18,18))#p_mat_m2_n1_m1_n2
		p_n1_m1_n2 = np.sum(p_m2_n1_m1_n2, axis = 0)
		p_m2_m1 = np.sum(np.sum(p_m2_n1_m1_n2, axis = 1), axis = 2)
		p_n1 = np.sum(np.sum(p_n1_m1_n2, axis = 1), axis = 1)
		p_n2 = np.sum(np.sum(p_n1_m1_n2,axis = 0), axis = 0)
		p_m1 = np.sum(p_m2_m1, axis = 1)
		p_m2 = np.sum(p_m2_m1, axis = 0)
		plt.plot(p_m2,label='q='+str(q)+',g='+str(g) )
plt.legend()
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