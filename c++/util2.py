import numpy
from numpy import *

import scipy.stats as stat
from scipy.stats import beta, binom
from scipy.special import gammaln
from math import exp, log

import csv

from data import Datum

from tssb import *

import traceback

def log_factorial(n):
	return gammaln(n + 1)

def log_bin_coeff(n, k):
	return gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1)

def log_binomial_likelihood(x, n, mu):	
	return x * log(mu) + (n - x) * log(1 - mu)

def log_beta(a, b):
	return gammaln(a) + gammaln(b) - gammaln(a + b)

def logsumexp(X, axis=None):
    maxes = numpy.max(X, axis=axis)
    return log(numpy.sum(exp(X - maxes), axis=axis)) + maxes

def load_data(fname):
	reader = csv.DictReader(open(fname), delimiter='\t')
	data = []
	id=0
	for row in reader:
		name = row['gene'] 
		a = [int(x) for x in row['a'].split(',')]       
		d = [int(x) for x in row['d'].split(',')]
        
		mu_r = [float(x) for x in row['mu_r'].split(',')]
		delta_r = [float(x) for x in row['delta_r'].split(',')]
        
		mu_v = [float(x) for x in row['mu_v'].split(',')]
		delta_v = [float(x) for x in row['delta_v'].split(',')]
        
		data.append(Datum(name, id, a, d, mu_r, mu_v, delta_r, delta_v))
		id+=1
	
	return data

def slice_sample2(init_x, logprob, sigma=.1, step_out=True, max_steps_out=1000,bounds=[0,1]):
	def dir_logprob(z): return logprob(z)
	r=numpy.random.rand()
	upper = min(bounds[1],init_x+ (1-r)*sigma)
	lower = max(bounds[0],init_x - r*sigma)
	llh_s = numpy.log(numpy.random.rand()) + dir_logprob(init_x)  
	l_steps_out = 0
	u_steps_out = 0
	if step_out:
		while dir_logprob(lower) > llh_s and l_steps_out < max_steps_out:			
			l_steps_out += 1
			lower -= sigma; lower=check_bounds(lower,u=bounds[1]);#shankar
		while dir_logprob(upper) > llh_s and u_steps_out < max_steps_out:
			u_steps_out += 1
			upper += sigma; upper=check_bounds(upper,u=bounds[1]);#shankar
	while True:
		new_z     = check_bounds((upper - lower)*numpy.random.rand()+lower) 
		new_llh   = dir_logprob(new_z)
		if numpy.isnan(new_llh):
			raise Exception("Slice sampler got a NaN")
		if new_llh > llh_s:
			break
		else:
			if new_z < init_x:
				lower = new_z
			elif new_z > init_x:
				upper = new_z
			else:
				raise Exception("Slice sampler shrank to zero!")
	return new_z  

def check_bounds(p,l=0.0001,u=.9999):
	if p < l: p=l
	if p > u: p=u
	return p