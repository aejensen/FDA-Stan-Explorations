data {
  int<lower = 0> N;  //Number of samples
  int<lower = 0> D;  //Number of observation points
  int<lower = 0> K;  //Number of components for reconstruction
  
  matrix[N, D] X;    //Functional data
}

transformed data {
  vector[D] mu;           //cross-sectional average  
  matrix[N, D] X_center;  //centered data
  
  for(d in 1:D) {
    mu[d] = mean(X[,d]); 
  }
  
  X_center = X - rep_matrix(mu, N)';
}

generated quantities {
  vector[D] lambda;        //eigenvalues
  matrix[D, D] psi;        //eigenvectors
  matrix[N, D] zeta;       //scores
  matrix[N, D] Xhat;       //reconstructed data
  
  vector[D] average = mu;  //cross-sectional average  
  
  {
    int order[D];
    matrix[D, D] covMat;
    
    covMat = (X_center' * X_center) ./ (N - 1);
    
    //Get eigendecomposition of covariance matrix
    lambda = eigenvalues_sym(covMat);
    psi = eigenvectors_sym(covMat);
    
    //Order decomposition in descending order of variance explained
    order = sort_indices_desc(lambda);
    lambda = lambda[order]; 
    psi = psi[, order];
    
    //Calculate scores
    zeta = X_center * psi;
    
    //Reconstruct data with 1:K PCs
    Xhat = rep_matrix(mu, N)' + zeta[, 1:K] * psi[, 1:K]';
  }
}
