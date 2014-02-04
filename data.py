from numpy import *
import scipy.stats as stat
from scipy.special import gammaln
import util2 as u

class Datum(object):
	def __init__(self, name,a, d, mu_r, mu_v, delta_r, delta_v):
		self.name = name
		self.a=a
		self.d=d
			
		self.mu_r=mu_r
		self.mu_v=mu_v
		
		self.delta_r = delta_r
		self.delta_v = delta_v
		
		self.log_pi_r = self._get_log_mix_wts(delta_r)
		self.log_pi_v = self._get_log_mix_wts(delta_v)
		
		self._log_bin_norm_const = [u.log_bin_coeff(self.d[tp], self.a[tp]) for tp in arange(len(self.a))]		
		
	def _get_log_mix_wts(self,delta):
		log_den = gammaln(sum(delta)+1)
		log_mix_wts = []
		for i, d_i in enumerate(delta):
			log_num = gammaln(d_i+1)
			for j, d_j in enumerate(delta):
				if i!=j: log_num += gammaln(d_j)
			log_mix_wts.append(log_num-log_den)
		return log_mix_wts
		
	# for multiple samples
	def _log_likelihood(self, phi):
		ntps = len(phi) # multi sample
		return sum([self.__log_likelihood__(phi[tp],tp) for tp in arange(ntps)])

	def __log_likelihood__(self, phi,tp):		
		ll = []        
		for mu_r, log_pi_r in zip(self.mu_r, self.log_pi_r):
			for mu_v, log_pi_v in zip(self.mu_v, self.log_pi_v):
				ll.append(self.__log_complete_likelihood__(phi, mu_r, mu_v,tp) + log_pi_r + log_pi_v)
		return u.logsumexp(ll)
    
	# for multiple samples
	def _log_complete_likelihood(self, phi, mu_r, mu_v):
		ntps = len(self.a)
		return sum([self.__log_complete_likelihood__(phi, mu_r, mu_v,tp) for tp in arange(ntps)])
		
	def __log_complete_likelihood__(self, phi, mu_r, mu_v,tp):
		mu = (1 - phi) * mu_r + phi * mu_v		
		return u.log_binomial_likelihood(self.a[tp], self.d[tp], mu) +  self._log_bin_norm_const[tp]


