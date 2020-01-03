mu <- function(t) 1*(t + sin(4*pi*t^2)^2)

psi1 <- function(t) sqrt(2) * sin(2*pi*t)

psi2 <- function(t) sqrt(2) * cos(2*pi*t)

#Check orthonormality
round(c(integrate(function(x) psi1(x)^2, 0, 1)$value,
        integrate(function(x) psi2(x)^2, 0, 1)$value,
        integrate(function(x) psi1(x) * psi2(x), 0, 1)$value), 6)

simData <- function(N, t, withMu, sigma) {
  simFunc <- function(t, withMu, sigma) {
    zeta1 <- rnorm(1, 0, 1)
    zeta2 <- rnorm(1, 0, 0.5)
    
    out <- zeta1 * psi1(t) + zeta2 * psi2(t)
    
    if(withMu) {
      out <- out + mu(t)
    }
    outNoise <- out + rnorm(length(t), 0, sigma)
    list(X = out, Xobs = outNoise)
  }
  
  d <- lapply(1:N, function(i) simFunc(t, withMu, sigma))
  dTrue <- t(sapply(d, function(q) q$X))
  dObs <- t(sapply(d, function(q) q$Xobs))
  list(true = dTrue, obs = dObs)
}

seCov <- function(s, t, sigma, l) {
  sigma^2 * exp(-(s-t)^2 / (2*l^2))
}

simData_long <- function(N, t, sigma) {
  zeta1 <- as.vector(rmvnorm(1, rep(0, N), outer(1:N, 1:N, seCov, sigma=8, l = N/5)))
  zeta2 <- as.vector(rmvnorm(1, rep(0, N), outer(1:N, 1:N, seCov, sigma=4, l = N/5)))
  zeta1 <- zeta1 - mean(zeta1)
  zeta2 <- zeta2 - mean(zeta2)

  dY <- lapply(1:N, function(i) {
    mu(tSeq) + zeta1[i] * psi1(tSeq) + zeta2[i] * psi2(tSeq)  
  })

  dY_noise <- lapply(1:N, function(i) {
    dY[[i]] + rnorm(length(t), 0, sigma)
  })
  
  dTrue <- t(sapply(dY, function(q) q))
  dObs <- t(sapply(dY_noise, function(q) q))  
  
  list(true = dTrue, obs = dObs, zeta1 = zeta1, zeta2 = zeta2)
}
