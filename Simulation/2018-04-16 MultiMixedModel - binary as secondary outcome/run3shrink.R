setwd("C:/Users/eandrinopoulou/Desktop")

library(JMbayes)
library(lme4)
library(rjags)

M = 10

beta1 <- matrix(, M, 5)
beta2 <- matrix(, M, 3)
beta3 <- matrix(, M, 3)
tauu <- numeric()
inv_D <- matrix(, M, 36)

for (i in 1:M){
  print(i)
  try(source("simulate3shrink.R"))
  try(source("analysis3shrink.R"))
  
  beta1[i, ] <- betas1
  beta2[i, ] <- betas2
  beta3[i, ] <- betas3
  tauu[i] <- tau
  inv_D[i, ] <- invD
}
