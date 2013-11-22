# Stochastic Simulation Algorithm input file
# Species n regulates itself with creation rate g(n) and per capita deathrate 1. Species m is regulated by species n with rate q(n) and per capita deathrate rho

# Reactions
Diffn1:
	n1 > n2
	D*n1
	
Deathn1:
	n1 > $pool
	n1
	
Birthn1:
	n1 > {2.0}n1
	(qminus*n0**nu + qplus*n1**nu)/(n1**nu + n0**nu)
	
Diffn2:
	n2 > n1
	D*n2
	
Deathn2:
	n2 > $pool
	n2
	
Birthn2:
	n2 > {2.0}n2
	(qminus*n0**nu + qplus*n2**nu)/(n2**nu + n0**nu)
	
Diffm1:
	m1 > m2
	D*m1
	
Deathm1:
	m1 > $pool
	m1
	
Birthm1:
	m1 > {2.0}m1
	(qminus*n0**nu + qplus*n1**nu)/(n1**nu + n0**nu)

Diffm2:
	m2 > m1
	D*m2
	
Deathm2:
	m2 > $pool
	m2
	
Birthm2:
	m2 > {2.0}m2
	(qminus*n0**nu + qplus*n2**nu)/(n2**nu + n0**nu)

# Fixed species

# Variable species
n1 = 0.0
n2 = 0.0
m1 = 0.0
m2 = 0.0

# Parameters
D = 1.0
rho = 1.0
qplus = 12.0
qminus = 1.0
n0 = 4.0
nu = 4.0
