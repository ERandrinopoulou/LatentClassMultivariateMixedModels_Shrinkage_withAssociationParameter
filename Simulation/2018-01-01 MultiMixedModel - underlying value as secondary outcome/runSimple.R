setwd("G:/2017-01 ROI/2017-01-03/Simulation")

library(JMbayes)
library(lme4)
library(rjags)

M = 2

beta1 <- matrix(, M, 3)
beta2 <- matrix(, M, 3)
tauu <- numeric()
inv_D <- matrix(, M, 16)
#alphaa <- mnumeric()

for (i in 1:M){
  print(i)
  source("simulateSimple.R")
  source("analysisSimple.R")
  
  beta1[i, ] <- betas1
  beta2[i, ] <- betas2
  tauu[i] < tau
  inv_D[i, ] <- invD
  #alphaa[i] <- alpha
}
