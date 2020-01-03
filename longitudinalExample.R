rm(list=ls())
library(mvtnorm)
library(fdapace)

source("util.R")

N <- 200
D <- 50
K <- 2
tSeq <- seq(0, 1, length.out=D)

set.seed(123456)
d <- simData_long(N, tSeq, 0.5)

matplot(tSeq, t(d$true), type="l", lty=1)
matplot(tSeq, t(d$obs), type="l", lty=1)

##########################################################
# Ordinary FPCA
##########################################################
m_fpca <- FPCA(lapply(1:N, function(i) d$true[i,]),
               lapply(1:N, function(i) tSeq))

plot(m_fpca)

plot(tSeq, -m_fpca$phi[,1])
lines(tSeq, psi1(tSeq))

plot(tSeq, -m_fpca$phi[,2])
lines(tSeq, psi2(tSeq))

plot(1:N, -m_fpca$xiEst[,1])
lines(1:N, d$zeta1)
plot(1:N, -m_fpca$xiEst[,2])
lines(1:N, d$zeta2)

