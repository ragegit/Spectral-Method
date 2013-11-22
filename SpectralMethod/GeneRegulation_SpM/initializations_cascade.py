import numpy as np
import math

# length of cascade
L = 3

# cutoff for particle number and new basis elements:
n1max = 40
n2max = 40
n3max = 40
j1max = 40
j2max = 40
j3max = 40

# system rates:
rho2 = 1.
rho3 = 1.

# other constants:
#g = 8. #shouldn't be a constant!
qplus = 12.
qminus = 1.
n0 = 4. #in paper
nu = 4.
# constant input creation rate
g = 8.