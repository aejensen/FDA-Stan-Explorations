/*
Model:
     mu | beta_mu = basis_mu * beta_mu
       b | beta_b = basis_b * beta_b
            eta_i = b * c_i
X_i(t) | eta_i(t) ~ N(mu(t) + eta_i(t), sigma)
*/

functions {
  real trapz(vector t, vector y) {
    /* Numerical integration of y(t) using the trapezoidal rule */
    /* Should ideally first check that t is sorted */
    real trapzsum = 0;
    int len = num_elements(t) - 1;
    for (ind in 1:len) {
      trapzsum += 0.5 * (t[ind + 1] - t[ind]) * (y[ind] + y[ind + 1]); 
    } 
    return trapzsum;
  }
  
  real norm2(vector t, vector y) {
    return sqrt(trapz(t, y .* y));
  }
}

data {
  int<lower = 1> N;              //Number of curves
  int<lower = 1> D;              //Number of sampling points
  int<lower = 1> K;              //Number of subject-level components
  int<lower = 1> P;              //Number of evaluation points
  
  vector[P] tEval;

  int<lower = 1> P_mu;           //Size of spline basis for population mean
  matrix[D, P_mu] basis_mu;      //Spline basis for mean function
	matrix[P, P_mu] basisEval_mu;  //Spline basis for evaluation of mean function
	
  int<lower = 1> P_b;            //Size of subject-level spline basis
  matrix[D, P_b] basis_b;        //Subject-level spline basis
  matrix[P, P_b] basisEval_b;    //Evaluation subjec-level spline basis
	
  vector[D] X[N];                //Functional data on regular grid
}

parameters {
  vector[P_mu] beta_mu;          //Population mean spline coefficients
  matrix[K, P_b] beta_b;         //Subject-level spline coefficients
  vector[K] c[N];                //Subject-level random effects
  real<lower = 0> sigma;         //Residual standard deviation
}

transformed parameters {
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
  {
    vector[D] mu = basis_mu * beta_mu;
    matrix[D, K] b = basis_b * beta_b';
    for (i in 1:N) {
      X[i] ~ normal(mu + b * c[i], sigma);
    }
  }
}

generated quantities {
  vector[P] mu;
  matrix[P, K] phi;
  vector[K] xi[N];
  vector[P] Xhat[N];
  
  mu = basisEval_mu * beta_mu;
  
  {
    matrix[P, K] bEval = basisEval_b * beta_b';
    phi = mdivide_right_tri_low(bEval, cholesky_decompose(bEval' * bEval));
    for(k in 1:K) {
      phi[, k] = phi[, k] / norm2(tEval, phi[, k]);
    }
  
    for(i in 1:N) {
      for(k in 1:K) {
        xi[i][k] = trapz(tEval, (bEval * c[i]) .* phi[, k]);
      }
    }
  }
  
  for(i in 1:N) {
    Xhat[i] = mu + (phi * xi[i]);
  }
}
