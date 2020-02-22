data {
  int<lower=0> n; // H_num
  int<lower=0> s; // Number of species
  vector[n] y; // data points
  matrix [n,s] M; // not sure
}
parameters {
  vector[s] beta; // species effects
  real<lower=0> sigma;
}
model {
  vector[n] mu; // general mean? 
  
  //Priors
  beta ~ normal(0, 10);
  sigma ~ normal(0,10);
  
  //Likelihood
  mu = M * beta;
  y ~ normal(mu, sigma);
}


