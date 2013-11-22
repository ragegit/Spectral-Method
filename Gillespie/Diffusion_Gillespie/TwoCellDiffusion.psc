# Stochastic Simulation Algorithm input file
# Species n regulates itself with creation rate g(n) and per capita deathrate 1. Species m is regulated by species n with rate q(n) and per capita deathrate rho

# Reactions
R1:
	n1 > n2
	D1*n1

R2:
	n2 > n1
	D2*n2

# Fixed species

# Variable species
n1 = 50.0
n2 = 50.0

# Parameters
D1 = 1.0
D2 = 1.0
