from initializations_Gillespie import *
import numpy as np
# functions for the rates:


def hill(n,nu,n0):
	""" hill function """
	return np.float(n**nu)/(n**nu + n0**nu)

def q_step(n,n0):
	""" n dependend creation rate for output species """
	#return qminus + (qplus - qminus)*hill(n,nu,n0)
	if n<=n0:
		return qminus
	else:
		return qplus
	
def q_hill(n,n0):
	""" n dependend creation rate for output species """
	#return qminus + (qplus - qminus)*hill(n,nu,n0)
	return np.float(qminus*n0**nu + qplus*n**nu)/(n**nu + n0**nu)
	
def q_n0(n):
	""" auxiliary function for q with one argument"""
	return q_hill(n,n0)

def g_n0(n):
	""" auxiliary function for g with one argument. For this example species n regulates itself with the same rate as species m is regulated by n"""
	return q_hill(n,n0)


