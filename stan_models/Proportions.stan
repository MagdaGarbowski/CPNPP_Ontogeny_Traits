data{
  int N; // number of observations
  int K; // number of groups in x (H_num)
  int x[N]; // observations at each x
  vector[N] y; // data points - MATT - would vector y[N] be the same? 
  int J; // number of populations
  int <lower = 1, upper = J> p[N]; // population ID 
}

parameters{
  vector<lower = 0, upper = 1> [K] betas; // intercept and slope for H_num 
  real<lower = 0> sigma; // error sd 
  vector[J] u; // POP_ID intercepts
  real<lower = 0> sigma_u; // POP_ID sd
}

model{
  /* Same as below - just not vectorized
	 for(n in 1:N)
	     y[n] ~ normal(betas[x[n]], sigma);
  */
  u ~ normal (0, sigma_u); // random effect POP_ID
  y ~ normal(betas[x] + u[p], sigma);
}
