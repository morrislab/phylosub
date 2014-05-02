#include <iostream>
#include<fstream>
#include<string>
#include<time.h>
#include<math.h>
#include<cstdlib>
#include<cstring>
#include <sstream>
#include<map>

#include "mh.hpp"
#include "util.hpp"


using namespace std;


//  g++ -o mh.o  mh.cpp  util.cpp `gsl-config --cflags --libs`
// ./mh.o 5000 0.01 14 1 17 6
//https://www.gnu.org/software/gsl/manual/html_node/Shared-Libraries.html


int main(int argc, char* argv[]){
	
	// parse command line args
	struct config conf;
	conf.MH_ITR=atoi(argv[1]);//5000
	conf.MH_STD=atof(argv[2]);//0.01
	conf.NDATA=atoi(argv[3]);//14; // no. of data points
	conf.NDELTA=atoi(argv[4]); //2, no. of genotypes
	conf.NNODES=atoi(argv[5]); //10, no. of nodes in the tree
	conf.TREE_HEIGHT=atoi(argv[6]);//5 
	
	char* FNAME_DATA = argv[7];
	
	struct datum data[conf.NDATA];
	load_data(FNAME_DATA,data,conf);
	for(int i=0;i<conf.NDATA;i++)
		data[i].log_pi_r=0;
	
	struct node nodes[conf.NNODES];
	load_tree("c_tree.txt",nodes);	 
	
	//start MH loop	
	mh_loop(nodes,data,conf);	
	
	// write updated params to disk
	write_params("c_params.txt",nodes,conf);	
	
	return 0;	
}


void mh_loop(struct node nodes[],struct datum data[],struct config conf){
	gsl_rng *rand = gsl_rng_alloc(gsl_rng_mt19937);
	for (int itr=0;itr<conf.MH_ITR;itr++){
		sample_cons_params(nodes,conf,rand);		
		
		double pi[conf.NNODES],pi_new[conf.NNODES];
		get_pi(nodes,pi_new,conf,0);
		get_pi(nodes,pi,conf,1);		
		
		double a = param_post(nodes,data,0,conf)-param_post(nodes,data,1,conf);			
		
		// apply the dirichlet correction terms
		double theta[conf.NNODES];// dirichlet params
		for(int i=0;i<conf.NNODES;i++)
			theta[i]=conf.MH_STD*pi_new[i];
		a += gsl_ran_dirichlet_lnpdf(conf.NNODES,theta,pi);
			
		for(int i=0;i<conf.NNODES;i++)
			theta[i]=conf.MH_STD*pi[i];
		a -= gsl_ran_dirichlet_lnpdf(conf.NNODES,theta,pi_new);
				
		double r = gsl_rng_uniform_pos(rand);		
		if (log(r)<a){
			update_params(nodes,conf);			
		}
	}
	gsl_rng_free(rand);
}


void sample_cons_params(struct node nodes[], struct config conf, gsl_rng *rand){

	map <int, int> node_id_map;
	for(int i=0;i<conf.NNODES;i++)
		node_id_map[nodes[i].id]=i;
		
	int NNODES=conf.NNODES;
	double pi[NNODES];
	for(int i=0;i<NNODES;i++)
		pi[i]=nodes[i].pi;
	
	// randomly sample from a dirichlet
	double pi_new[NNODES],alpha[NNODES];
	for(int i=0;i<NNODES;i++)
		alpha[i]=conf.MH_STD*pi[i]+1;
	dirichlet_sample(NNODES,alpha,pi_new,rand);
	
	// update the nodes pi1 (new pi)
	for(int i=0;i<NNODES;i++)
		nodes[i].pi1=pi_new[i];		
	
	// update the nodes param1 (new param)
	for(int i=0;i<NNODES;i++){		
		double param = nodes[i].pi1;			
		for(int c=0;c<nodes[i].nchild;c++){
			param+=nodes[node_id_map[nodes[i].cids.at(c)]].param1;
		}
		nodes[i].param1=param;
	}
}


double param_post(struct node nodes[], struct datum data[], int old,struct config conf){	
	double llh = 0.0;
	for(int i=0;i<conf.NNODES;i++){
		double p=0;
		if(old==0)
			p=nodes[i].param1;
		else
			p=nodes[i].param;
		for(int j=0;j<nodes[i].ndata;j++){
			llh+=data[nodes[i].dids.at(j)].log_ll(p);
		}
	}
	return llh;	
}

void update_params(struct node nodes[],struct config conf){	
	for(int i=0;i<conf.NNODES;i++){
		nodes[i].param=nodes[i].param1;
		nodes[i].pi=nodes[i].pi1;
	}
}

void get_pi(struct node nodes[], double pi[], struct config conf, int old){
	for(int i=0;i<conf.NNODES;i++){
		if (old==0)
			pi[i]=nodes[i].pi1;
		else
			pi[i]=nodes[i].pi;
	}
}


void load_data(char fname[],struct datum *data, struct config conf){
	string line,token,token1;
	ifstream dfile (fname);
	int ctr=0,id=-1;
	while (getline (dfile,line,'\n')){
		if (id==-1){id+=1;continue;}
		istringstream iss(line);
		ctr=0;
		while(getline(iss,token,'\t')){
			if(ctr==0){
				data->id=id;
			}
			else if(ctr==2){
				data->a=atoi(token.c_str());
			}
			else if(ctr==3){
				data->d=atof(token.c_str());
				data->log_bin_norm_const=log_bin_coeff(data->d,data->a);
			}
			else if(ctr==4){
				data->mu_r=atof(token.c_str());
			}
			else if(ctr==5){
				data->delta_r=atoi(token.c_str());
			}
			else if(ctr==6){
				istringstream iss(token);
				for(int i=0;i<conf.NDELTA;i++){
					getline(iss,token1,',');
					data->mu_v.push_back(atof(token1.c_str()));
				}
			}
			else if(ctr==7){
				istringstream iss(token);
				for(int i=0;i<conf.NDELTA;i++){
					getline(iss,token1,',');
					data->delta_v.push_back(atoi(token1.c_str()));
				}
				data->set_log_mix_wts(data->delta_v);
			}
			ctr+=1;			
		}		
		data++;
		id+=1;
	}	
	dfile.close();	
}


void load_tree(char fname[], struct node *nodes){
	string line,token,token1;
	ifstream dfile (fname);
	int ctr=0;
	while (getline (dfile,line,'\n')){
		istringstream iss(line);
		ctr=0;
		while(getline(iss,token,'\t')){
			if(ctr==0){
				nodes->id=atoi(token.c_str());
			}
			else if(ctr==1){
				nodes->param=atof(token.c_str());
			}
			else if(ctr==2){
				nodes->pi=atof(token.c_str());
			}
			else if(ctr==5){
				nodes->ndata=atoi(token.c_str());
			}
			else if(ctr==6){
				istringstream iss(token);
				for(int i=0;i<nodes->ndata;i++){
					getline(iss,token1,',');
					nodes->dids.push_back(atoi(token1.c_str()));
				}
			}
			else if(ctr==3){
				nodes->nchild=atoi(token.c_str());		
			}			
			else if(ctr==4){
				istringstream iss(token);
				for(int i=0;i<nodes->nchild;i++){
					getline(iss,token1,',');
					nodes->cids.push_back(atoi(token1.c_str()));
				}
			}
			else if(ctr==7){
				nodes->ht=atoi(token.c_str());		
			}			
			ctr+=1;
		}
		nodes++;
	}	
	dfile.close();	
}

void write_params(char fname[], struct node *nodes, struct config conf){	
	ofstream dfile;
	dfile.open(fname);
	for(int i=0;i<conf.NNODES;i++){		
		dfile<<nodes[i].id<<'\t'<<nodes[i].param<<'\t'<<nodes[i].pi<<'\n';
	}	
	dfile.close();	
}