

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
  int<lower = 1, upper = n_pop> sp_Hnum[N]; // populations
  vector[N] y; // data points
}
transformed data{
  // For HS prior
  real<lower = 0> scale_global = 1; // scale for half-t for tau
  real<lower = 1> nu_global = 1; //df for half-t prior for tau
  real<lower = 1> nu_local = 1; //df for half-t for lambda
}

parameters {
  real alpha; 
  vector[n_Hnum - 1] beta_Hnum_raw;
  vector[n_sp - 1] beta_sp_raw;
  vector[n_pop - 1] beta_pop_raw; 
  vector[n_pop - 1] beta_sp_Hnum_raw; 
  real<lower=0> sigma;
  // For HS prior
  real<lower = 0> r1_global;
  real<lower = 0> r2_global;
  vector<lower = 0>[4] r1_local;
  vector<lower = 0>[4] r2_local;
}

transformed parameters{
  vector[n_Hnum] beta_Hnum;
  vector[n_sp]   beta_sp;
  vector[n_pop] beta_pop; 
  vector[n_pop] beta_sp_Hnum;  
  // For HS prior
  real<lower = 0> tau = r1_global * sqrt(r2_global);
  vector[4] lambda = r1_local .* sqrt(r2_local);
  vector<lower = 0>[4] raneff_sigmas = tau * lambda;

  beta_Hnum = raw_to_beta(beta_Hnum_raw, raneff_sigmas[1]);      
  beta_sp = raw_to_beta(beta_sp_raw, raneff_sigmas[2]);          
  beta_pop = raw_to_beta(beta_pop_raw, raneff_sigmas[3]);         
  beta_sp_Hnum = raw_to_beta(beta_sp_Hnum_raw, raneff_sigmas[4]); 
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

  // For HS prior
  r1_local ~ normal(0, 1);
  r2_local ~ inv_gamma(0.5 * nu_local, 0.5 * nu_local);
  r1_global ~ normal(0, scale_global * sigma);
  r2_global ~ inv_gamma(0.5 * nu_global, 0.5 * nu_global);
  
  sigma ~ normal(0,10);
  
  //Likelihood
  y ~ normal(mu, sigma);
}


