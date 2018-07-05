setwd("C:/Users/eandrinopoulou/Desktop")

library(JMbayes)
library(lme4)
library(rjags)

M = 10

beta1 <- matrix(, M, 4)
beta2 <- matrix(, M, 3)
tauu <- numeric()
inv_D <- matrix(, M, 16)

for (i in 1:M){
  print(i)
  try(source("pbc_simulate.R"))
  try(source("pbc_analysis.R"))
  
  beta1[i, ] <- betas1
  beta2[i, ] <- betas2
  tauu[i] <- tau
  inv_D[i, ] <- invD
}
