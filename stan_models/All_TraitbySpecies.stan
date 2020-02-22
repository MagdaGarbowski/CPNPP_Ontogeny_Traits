// Means mod: y = alpha (species) + beta (H_num)  + alpha * beta + error //
// Effects mod: y = mu + alpha (species) + beta (H_num) + alpha * beta + error//
data{
  int N; // number of observations
  int K; // number of H_num levels
  int S; // number of Species
  int<lower = 0, upper = K> x[N]; // H_num observations ID
  int<lower = 1, upper = S> s[N]; // Species observations ID
  vector[N] y; // data points
  int J; // number of populations
  int<lower = 1, upper = J> p[N]; // population ID 
  
}

parameters{
  vector<lower = 0> [K] betas; // intercepts and slopes for H_num 
  vector<lower = 0> [S] alphas; // intercepts and slopes for Species
  real<lower = 0> sigma; // error sd 
  vector[J] u_raw; // POP_ID intercepts
  real<lower = 0> sigma_u; // POP_ID sd
}

transformed parameters{
  /* 
	 u_raw ~ normal(0,1)
	 u = u_raw * sigma_u 
	 
	 is the same as u ~ normal(0, sigma_u) 
	 but sometimes samples more efficently
  */
  
  vector[J] u = u_raw * sigma_u;
}

model{
  /* Questions: (1)
  (1) How to code interaction? 
  (2) How to add random effect of POP_ID? 
  */
  // weak priors on sigmas
  sigma ~ normal(0, 5);
  sigma_u ~ normal(0, 5);
  u_raw ~ normal(0, 1); // random effect POP_ID - non-centered parameterization
  
  	 for(n in 1:N)
	     y[n] ~ normal(betas[x[n]] + alphas[s[n]] + betas[x[n]] * alphas[s[n]], sigma);
}



