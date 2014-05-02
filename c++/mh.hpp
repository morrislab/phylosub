
#include<vector>
#include<math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


#include "util.hpp"


void sample_cons_params(struct node nodes[],struct config conf,gsl_rng *rand);
double param_post(struct node nodes[], struct datum data[], int old,struct config conf);
void update_params(struct node nodes[], struct config conf);
void get_pi(struct node nodes[], double pi[], struct config conf, int old);

void load_data(char fname[], struct datum data[], struct config conf);
void load_tree(char fname[], struct node nodes[]);
void write_params(char fname[], struct node nodes[], struct config conf);

void mh_loop(struct node nodes[], struct datum data[], struct config conf);

struct config{
	int MH_ITR;
	float MH_STD;
	
	int NDATA; // no. of data points
	int NDELTA; // no. of genotypes
	int NNODES; // no. of nodes in the tree
	int TREE_HEIGHT; 	
	
};

struct datum{
	
	int id;	
	int a,d;
	
	double mu_r;
	vector<double> mu_v; // double mu_v[NDELTA];
	int delta_r;
	vector<int> delta_v; // int delta_v[NDELTA] 
	
	double log_pi_r;// this is just 0 for scalar. 
	vector<double> log_pi_v; //double log_pi_v[NDELTA]; //set_log_mix_wts()
	
	double log_bin_norm_const;//log_bin_coeff(d,a);	
	
	void set_log_mix_wts(vector<int> delta){
		int sd=0;
		int NDELTA = delta.size();
		for(int i=0;i<NDELTA;i++)
			sd+=delta[i];
		double log_den = lgamma(sd+1);
		
		for(int i=0;i<NDELTA;i++){
			double log_num = lgamma(delta[i]+1);
			for(int j=0;j<NDELTA;j++)
				if (i!=j)
					log_num += lgamma(delta[j]);			
			log_pi_v.push_back(log_num-log_den);
		}
	}
	
	
	double log_ll(double phi){
		int NDELTA = delta_v.size();
		double ll[NDELTA];
		for(int i=0;i<NDELTA;i++){
			ll[i] = log_complete_ll(phi,mu_r,mu_v[i]) + log_pi_r + log_pi_v[i];
		}
		return logsumexp(ll,NDELTA);
	}
	
	double log_complete_ll(double phi, double mu_r, double mu_v){		
		double mu = (1 - phi) * mu_r + phi * mu_v;	
		return log_binomial_likelihood(a, d, mu) + log_bin_norm_const;	
	}
	
};


struct node{
	int id;
	double param,pi;
	double param1,pi1; // dummy	
	int ndata;
	vector<int> dids;	
	int nchild;
	vector<int> cids; // children ids	
	int ht;	
};