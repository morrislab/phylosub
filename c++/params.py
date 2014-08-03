
###### code to sample from the paramater posterior p(\phi | data) ########

import numpy
from numpy import *
from data import Datum

from tssb import *

from util import dirichletpdfln
from numpy.random import dirichlet

from subprocess import call

def metropolis(tssb,iters=1000,std=[0.01],burnin=0,ntps=1,fin=''):

	wts, nodes = tssb.get_mixture()

	for tp in arange(ntps):
		sample_cons_params(tssb,tp)
		update_params(tssb,tp)
	
	tp=0
	
	set_node_height(tssb)	
	write_tree(tssb)
		
	MH_ITR = str(iters)
	MH_STD = str(std[0])
	NDATA = str(tssb.num_data)
	NDELTA = str(len(tssb.data[0].delta_v))
	NNODES = str(len(nodes))
	TREE_HEIGHT = str(max([node.ht for node in nodes])+1)
	NTPS = str(ntps)
	
	FNAME_DATA = fin
	
	#call(['./mh.o', '-d', './vw/'+fname, '-c', '--passes', '25', '--l1', '1e-6', '--l2', '0', '--readable_model', './weights/'+fname])
	call(['./mh.o', MH_ITR, MH_STD, NDATA, NDELTA, NNODES, TREE_HEIGHT,FNAME_DATA, NTPS])
	
	update_tree_params(tssb)	
	
	return 0


def write_tree(tssb):
	fh=open('c_tree.txt','w')
	wts,nodes=tssb.get_mixture()
	
	def descend(root):		
		for child in root.children():			
			descend(child)
		
		# write data#
		cids=''
		for child in root.children():cids+=str(child.id)+','
		cids=cids.strip(',')
		if cids=='': cids=str(-1)
		
		dids=''
		for dat in root.get_data():dids+=str(dat.id)+','
		dids=dids.strip(',')
		if dids=='': dids=str(-1)
		
		line = str(root.id) + '\t' + list_to_string(root.params) + '\t' + list_to_string(root.pi) + '\t' + str(len(root.children())) + '\t'  + cids + '\t' + str(len(root.get_data())) + '\t' + dids + '\t' +  str(root.ht)
		fh.write(line)
		fh.write('\n')
		fh.flush()
		###############
	
	descend(tssb.root['node'])
	fh.flush()
	fh.close()

def list_to_string(p):
    o=''
    for pp in p:o+=str(pp)+','
    return o.strip(',')

def update_tree_params(tssb):
	wts, nodes = tssb.get_mixture()
	ndict = dict()
	for node in nodes: ndict[node.id]=node
	
	fh=open('c_params.txt')
	params=[line.split() for line in fh.readlines()]
	fh.close()
    
	for p in params:
		ndict[int(p[0])].params = string_to_list(p[1])
		ndict[int(p[0])].pi = string_to_list(p[2])

def string_to_list(p):
    p=p.strip(',')
    return array([float(pp) for pp in p.split(',')])
    
def set_node_height(tssb):
	tssb.root['node'].ht=0
	def descend(root,ht):
		for child in root.children():
			child.ht=ht
			descend(child,ht+1)
	descend(tssb.root['node'],1)
	
		
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
	
	
def update_params(tssb,tp):
	def descend(root,tp):			
		for child in root.children():
			descend(child,tp)	
		root.params[tp] = root.params1[tp]
		root.pi[tp] = root.pi1[tp]
	descend(tssb.root['node'],tp)