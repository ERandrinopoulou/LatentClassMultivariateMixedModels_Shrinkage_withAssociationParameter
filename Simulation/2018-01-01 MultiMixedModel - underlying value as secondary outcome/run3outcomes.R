setwd("C:/Users/015004/Desktop/2017-01 ROI/2017-01-03/Simulation")

library(JMbayes)
library(lme4)
library(rjags)

M = 200

beta1 <- matrix(NA, M, 3)
beta2 <- matrix(NA, M, 3)
beta3 <- matrix(NA, M, 3)
tauu <- rep(NA, M) 
inv_D <- matrix(NA, M, 36)
alpha1 <- rep(NA, M) 
alpha2 <- rep(NA, M) 

for (i in 5:M){
  print(i)
  try(source("simulate3outcomes.R"))
  try(source("analysis3outcomes.R"))
  
  beta1[i, ] <- betas1
  beta2[i, ] <- betas2
  beta3[i, ] <- betas3
  tauu[i] <- tau
  inv_D[i, ] <- invD
  alpha1[i] <- alphas1
  alpha2[i] <- alphas2
}