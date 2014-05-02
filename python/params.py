
###### code to sample from the paramater posterior p(\phi | data) ########

import numpy
from numpy import *
from data import Datum

from tssb import *

from util import dirichletpdfln
from numpy.random import dirichlet

def metropolis(tssb,iters=1000,std=[0.01],burnin=0,ntps=1):

	wts, nodes = tssb.get_mixture()
	pi = empty((ntps,len(wts)))

	for tp in arange(ntps):
		sample_cons_params(tssb,tp)
		update_params(tssb,tp)
		pi[tp] = getpi(tssb,tp)
		
	ctr=0
	pi_new  = empty(pi.shape)	
	for i in arange(-burnin,iters):
		cf=0
		for tp in arange(ntps): 
			pi_new[tp] = sample_cons_params1(tssb,pi[tp],tp,std)
		a=multi_param_posterior(tssb,1,ntps)-multi_param_posterior(tssb,0,ntps) + multi_correction_term(pi,pi_new,std[0],ntps) - multi_correction_term(pi_new,pi,std[0],ntps)
		if log(rand(1)) < a:
			for tp in arange(ntps): 
				update_params(tssb,tp)
				pi[tp] = pi_new[tp]
			ctr+=1		
	return ctr*1./iters


# tree-structured finite-dimensional stick breaking
def sample_cons_params(tssb,tp):
	def descend(root,tp):
	
		if root.parent() is None:
			root.params1[tp] = 1
			root.pi1[tp] = root.params1[tp]*rand(1) # break nu stick
		r = root.params1[tp]-root.pi1[tp] #mass assigned to children
		p = rand(len(root.children()));p=r*p*1./sum(p)
		index=0
		for child in root.children():			
			#print child.params1
			child.params1[tp] = p[index]# break psi sticks			
			child.pi1[tp] = child.params1[tp]*(rand(1)**(len(child.children())>0)) # break nu stick
			index+=1
		for child in root.children():			
			descend(child,tp)	

	descend(tssb.root['node'],tp)
	
			
# no stick breaking, updates pi with small perturbations
def sample_cons_params1(tssb,pi,tp,mh_std=[0.01]):	
	std = mh_std[randint(0,len(mh_std))]
	pi = dirichlet(std*pi+1) # for dirichlet proposal
	pi = pi+0.0001; pi=pi/sum(pi)
	pi=list(pi);pi.reverse()
		
	def descend(root,tp):		
		for child in root.children():			
			descend(child,tp)
		root.pi1[tp]=pi.pop();
		if len(root.children())==0:
			root.params1[tp] = root.pi1[tp]
		else:			
			root.params1[tp] = root.pi1[tp]+sum([child.params1[tp] for child in root.children()])			
	descend(tssb.root['node'],tp)
	return getpi(tssb,tp,1)

# MH correction terms for asymmetric proposal distribution
def multi_correction_term(pi1,pi2,std,ntps):
	return sum([dirichletpdfln(pi1[tp],std*pi2[tp]) for tp in arange(ntps)]) # for dirichlet proposal

def multi_param_posterior(tssb, new,ntps):
	return sum([param_posterior(tssb,new,tp) for tp in arange(ntps)])

def param_posterior(tssb, new,tp):
	def getposterior(root,new,tp):
		llh = 0
		data = root.get_data()
		if new: 
			p = root.params1[tp]
		else:
			p = root.params[tp]
		llh = sum([data[i].__log_likelihood__(p,tp) for i in arange(len(data))])	

		for child in root.children():					
			llh += getposterior(child,new,tp)
		return llh
	return getposterior(tssb.root['node'],new,tp)	

def update_params(tssb,tp):
	def descend(root,tp):			
		for child in root.children():
			descend(child,tp)	
		root.params[tp] = root.params1[tp]
		root.pi[tp] = root.pi1[tp]
	descend(tssb.root['node'],tp)

def getpi(tssb,tp,new=0):
	pi = []
	def descend(root,tp):
		for child in root.children():		
			descend(child,tp)
		if new: 
			p = root.pi1[tp]
		else:
			p = root.pi[tp]
		pi.append(p)
	descend(tssb.root['node'],tp)
	pi=array(pi);pi.shape=len(pi),	
	return pi
