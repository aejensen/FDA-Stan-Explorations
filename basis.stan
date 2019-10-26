/*
Model:
     mu | beta_mu = basis_mu * beta_mu
       b | beta_b = basis_b * beta_b
            eta_i = b * c_i
X_i(t) | eta_i(t) ~ N(mu(t) + eta_i(t), sigma)
*/

data {
  int<lower = 1> N;           //Number of curves
  int<lower = 1> D;           //Number of sampling points
  int<lower = 1> K;           //Number of subject-level components

  int<lower = 1> P_mu;        //Size of spline basis for population mean
  matrix[D, P_mu] basis_mu;   //Spline basis for population mean
	
  int<lower = 1> P_b;         //Size of subject-level spline basis
  matrix[D, P_b] basis_b;     //Subject-level spline basis
	
  vector[D] X[N];             //Functional data on regular grid
}

parameters {
  vector[P_mu] beta_mu;       //Population mean spline coefficients
  matrix[K, P_b] beta_b;      //Subject-level spline coefficients
  vector[K] c[N];             //Subject-level random effects
  real<lower = 0> sigma;      //Residual standard deviation
}

transformed parameters {
  vector[D] mu;               //Population average
  matrix[D, K] b;             //Basis functions for subject-level effects
  vector[D] eta[N];           //Subject-level functions  

  mu = basis_mu * beta_mu;
  
  b = basis_b * beta_b';
  for(i in 1:N) {
    eta[i] =  b * c[i];
  }
}

model {
  /* Spline coefficient priors */
  beta_mu ~ normal(0, 50);
  for(k in 1:K) {
    beta_b[k,] ~ normal(0, 50);
  }  

  /* Subject random effects prior */  
  for (i in 1:N) {
    c[i] ~ normal(0, 1);
  }

  /* Residual standard deviation prior */
  sigma ~ cauchy(0, 5);
  
  /* Likelihood */
  for (i in 1:N) {
    X[i] ~ normal(mu + eta[i], sigma);
  }
}

generated quantities {
  vector[D] Xhat[N]; // Fitted functional functions

  for(i in 1:N) {
    Xhat[i] = mu + eta[i];
  }
}

