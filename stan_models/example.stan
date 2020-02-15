data{
  real mu1; // mu  y 
  real sigma1; // sigma y 
  real mu2; // mu var
  real sigma2; // sigma var
  int N; // number of observations
  int K; // number of traits
  int P[N]; // Trait ID
  vector [N] y; //datapoints
}

parameters{
  vector[K] mu;
  vector<lower = 0>[K] sigma;
}

model{
  mu ~ normal(mu1, sigma1);
  sigma ~ normal(mu2, sigma2);
  y ~ normal(mu[P], sigma[P]);
}

generated quantities{
  vector[K] mu_out;
  vector[K] sigma_out;
  vector[K] cv; 
  for(i in 1:K){
  mu_out[i] = mu[i];
  sigma_out[i] = sigma[i];
	cv[i] = sigma_out[i]/ mu_out[i];
  }
  
}

