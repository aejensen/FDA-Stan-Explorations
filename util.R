mu <- function(t) 1*(t + sin(4*pi*t^2)^2)

psi1 <- function(t) sqrt(2) * sin(2*pi*t)

psi2 <- function(t) sqrt(2) * cos(2*pi*t)

#Check orthonormality
round(c(integrate(function(x) psi1(x)^2, 0, 1)$value,
        integrate(function(x) psi2(x)^2, 0, 1)$value,
        integrate(function(x) psi1(x) * psi2(x), 0, 1)$value), 6)

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

simData <- function(N, t, withMu, sigma) {
  d <- lapply(1:N, function(i) simFunc(t, withMu, sigma))
  dTrue <- t(sapply(d, function(q) q$X))
  dObs <- t(sapply(d, function(q) q$Xobs))
  list(true = dTrue, obs = dObs)
}
