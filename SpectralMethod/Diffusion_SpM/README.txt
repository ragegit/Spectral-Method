reproduce-Walczac
=================

Functionality:
Spectral method implemented in python for a two species cascade gene regulatory network.

Aim:
Gene regulation cascades at multiple sites are coupled via diffusion

Files:
- initializations.py
	contains all variables which do not change during one simulation
- functions.py
	contains functions to create:
		particle dependend rates
		input probabilities
		vectors and matrices for the final recursion equation of the spectral method
- main.py
	main code to run the simulation. After running main.py, the full probability density p(n,m) of a two species gene regulation cascade is calculated. Run main.py in ipython and plot the matrix p_mat in 3D to see p(n,m). The marginal distributions p(n), p(m) are called p_n and p_m respectively. Mean values are called n_mean and m_mean.
	
- spectral_plots.py:
	plots the full probability density p(n,m) and the marginal probability densities
