setwd("C:/Users/eandrinopoulou/Desktop")

library(JMbayes)
library(lme4)
library(rjags)

M = 1

beta11 <- matrix(, M, 3)
beta12 <- matrix(, M, 3)
beta21 <- matrix(, M, 3)
beta22 <- matrix(, M, 3)
beta31 <- matrix(, M, 3)
beta32 <- matrix(, M, 3)
alpha11 <- numeric()
alpha12 <- numeric()
alpha21 <- numeric()
alpha22 <- numeric()
tauu <- numeric()
inv_D1 <- matrix(, M, 36)
inv_D2 <- matrix(, M, 36)
prrr <- list()

for (i in 1:M){
  print(i)
  set.seed(i+1)
  
  ########################
  
  source("simulate3a_shrink.R")
  
  data_1 <- data
  data.id_1 <- data.id
  
  source("simulate3b_shrink.R")
  
  data_2 <- data
  data.id_2 <- data.id
  data_2$id <- data_2$id + tail(data.id_1$id, n=1)
  data.id_2$id <- data.id_2$id + tail(data.id_1$id, n=1)
  
  data <- rbind(data_1, data_2)
  data.id <- rbind(data.id_1, data.id_2)
  
  ########################
  
  try(source("analysis3shrink.R"))
  
  beta11[i, ] <- betas11
  beta12[i, ] <- betas12
  beta21[i, ] <- betas21
  beta22[i, ] <- betas22
  beta31[i, ] <- betas31
  beta32[i, ] <- betas32
  alpha11[i] <- alphas11 
  alpha12[i] <- alphas12 
  alpha21[i] <- alphas21 
  alpha22[i] <- alphas22 
  tauu[i] <- tau
  inv_D1[i, ] <- invD1
  inv_D2[i, ] <- invD2
  prrr[[i]] <- prr 
}


apply(beta11[c(1:14,18),],2,mean) ## simulate2a
apply(beta12[c(15:17),],2,mean) ## simulate2a

apply(beta11[c(15:17),],2,mean) ## simulate2b
apply(beta12[c(1:14,18),],2,mean) ## simulate2b


apply(beta21[c(1:14,18),],2,mean) ## simulate2a
apply(beta22[c(15:17),],2,mean) ## simulate2a

apply(beta21[c(15:17),],2,mean) ## simulate2b
apply(beta22[c(1:14,18),],2,mean) ## simulate2b


apply(beta31[c(1:14,18),],2,mean) ## simulate2a
apply(beta32[c(15:17),],2,mean) ## simulate2a

apply(beta31[c(15:17),],2,mean) ## simulate2b
apply(beta32[c(1:14,18),],2,mean) ## simulate2b


mean(alpha11[c(1:14,18)]) ## simulate2a
mean(alpha12[c(15:17)]) ## simulate2a

mean(alpha11[c(15:17)]) ## simulate2b
mean(alpha12[c(1:14,18)]) ## simulate2b



mean(alpha21[c(1:14,18)]) ## simulate2a
mean(alpha22[c(15:17)]) ## simulate2a

mean(alpha21[c(15:17)]) ## simulate2b
mean(alpha22[c(1:14,18)]) ## simulate2b
