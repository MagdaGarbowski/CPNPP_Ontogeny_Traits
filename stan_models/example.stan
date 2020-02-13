data{
  int N;
  int K; // number of groups in x
  int x[N];
  vector[N] y;
}

parameters{
  vector[K] betas;
  real<lower = 0> sigma;
}

model{
  
  /*
	Don't want to do this - we have 1 beta, and we assume that x is a
	continuous variable

	x is a group

	y ~ normal(beta * x, sigma);
  
  */

  /* Same as below - just not vectorized
	 for(n in 1:N)
	     y[n] ~ normal(betas[x[n]], sigma);
  */
  
  y ~ normal(betas[x], sigma);
}
