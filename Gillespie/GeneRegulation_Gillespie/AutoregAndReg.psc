# Stochastic Simulation Algorithm input file
# Species n regulates itself with creation rate g(n) and per capita deathrate 1. Species m is regulated by species n with rate q(n) and per capita deathrate rho

# Reactions
R1:
	n > {2} n
	(qminus*n0**nu + qplus*n**nu)/(n**nu + n0**nu)
	
R2:
	n > $pool
	n
    
R3:
	m > {2} m
	(qminus*n0**nu + qplus*n**nu)/(n**nu + n0**nu)
	
R4:
	m > $pool
	rho*m

# Fixed species
 
# Variable species
n = 10.0
m = 10.0

# Parameters
rho = 1.0
qplus = 12.0
qminus = 1.0
n0 = 4.0
nu = 4.0
g = 7.0
