setwd("C:/Users/015004/Desktop/2017-01 ROI/2017-01-03/Simulation")

library(JMbayes)
library(lme4)
library(rjags)

M = 2

beta1 <- matrix(NA, M, 3)
beta2 <- matrix(NA, M, 3)
tauu <- rep(NA, M) 
inv_D <- matrix(NA, M, 16)
alpha <- rep(NA, M) 


for (i in 1:M){
  print(i)
  try(source("simulate2outcomes2.R"))
  try(source("analysis2outcomes.R"))
  
  beta1[i, ] <- betas1
  beta2[i, ] <- betas2
  tauu[i] <- tau
  inv_D[i, ] <- invD
  alpha[i] <- alphas
}
