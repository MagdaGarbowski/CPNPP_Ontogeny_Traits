data{
  int N; // number of observations
  int K; // number of groups in x (H_num)
  int<lower = 0, upper = K> x[N]; // observations at each x
  vector[N] y; // data points
  int J; // number of populations
  int<lower = 1, upper = J> p[N]; // population ID 
}

parameters{
  vector<lower = 0> [K] betas; // intercept and slope for H_num 
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
  /* Same as below - just not vectorized
	 for(n in 1:N)
	     y[n] ~ normal(betas[x[n]], sigma);
  */
  // weak priors on sigmas
  sigma ~ normal(0, 5);
  sigma_u ~ normal(0, 5);

  u_raw ~ normal(0, 1); // random effect POP_ID - non-centered parameterization
  y ~ normal(betas[x] + u[p], sigma);
}
