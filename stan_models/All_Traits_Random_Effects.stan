

functions{

  // take any vector of raw_betas, and convert to beta
  vector raw_to_beta(vector raw_betas, real sigma){
	int N = dims(raw_betas)[1];
	vector[N + 1] betas = append_row(raw_betas, -sum(raw_betas)) * sigma;

	return betas;	
  }
	
}
data {
  int<lower = 1> N; // number of obs
  int<lower = 1> n_Hnum;
  int<lower = 1> n_sp;
  int<lower = 1> n_pop;
  int<lower = 1> n_sp_Hnum;
  int<lower = 1, upper = n_Hnum> Hnum[N]; // H_num
  int<lower = 1, upper = n_sp> sp[N]; // Species
  int<lower = 1, upper = n_pop> pop[N]; // populations
  int<lower = 1, upper = n_sp_Hnum> sp_Hnum[N]; // populations
  vector[N] y; // data points
}

parameters {
  real alpha; 
  vector[n_Hnum - 1] beta_Hnum_raw;
  vector[n_sp - 1] beta_sp_raw;
  vector[n_pop - 1] beta_pop_raw; 
  vector[n_sp_Hnum - 1] beta_sp_Hnum_raw; 
  vector<lower = 0>[4] raneff_sigmas;
  real<lower=0> sigma;
}

transformed parameters{
  vector[n_Hnum] beta_Hnum = raw_to_beta(beta_Hnum_raw, raneff_sigmas[1]);
  vector[n_sp] beta_sp = raw_to_beta(beta_sp_raw, raneff_sigmas[2]);
  vector[n_pop] beta_pop = raw_to_beta(beta_pop_raw, raneff_sigmas[3]); 
  vector[n_sp_Hnum] beta_sp_Hnum = raw_to_beta(beta_sp_Hnum_raw, raneff_sigmas[4]); 
}

model{
  vector[N] mu = alpha +
	beta_Hnum[Hnum] +
	beta_sp[sp] +
	beta_pop[pop] +
	beta_sp_Hnum[sp_Hnum]; // general mean? 

  //Priors
  beta_Hnum_raw ~ normal(0,1);
  beta_sp_raw ~ normal(0,1);
  beta_pop_raw ~ normal(0,1);
  beta_sp_Hnum_raw ~ normal(0,1);
  raneff_sigmas ~ cauchy(0, 0.1);
  sigma ~ normal(0,10);
  
  //Likelihood
  y ~ normal(mu, sigma);
}


