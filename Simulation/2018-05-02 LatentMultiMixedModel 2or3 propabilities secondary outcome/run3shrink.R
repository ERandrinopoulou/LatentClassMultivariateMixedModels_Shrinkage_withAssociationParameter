setwd("C:/Users/015004/Desktop/2018-04-18 LatentMiltiMixedModel 2or3 out binary")

library(JMbayes)
library(lme4)
library(rjags)

M = 20

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


apply(rbind(beta11[c(1:3,7,8,15,18,19),],beta12[c(4,5,6,10:13,16,20),]),2,mean) ## simulate3a
apply(rbind(beta12[c(1:3,7,8,15,18,19),],beta11[c(4,5,6,10:13,16,20),]),2,mean) ## simulate3b

apply(rbind(beta21[c(1:3,7,8,15,18,19),],beta22[c(4,5,6,10:13,16,20),]),2,mean) ## simulate3a
apply(rbind(beta22[c(1:3,7,8,15,18,19),],beta21[c(4,5,6,10:13,16,20),]),2,mean) ## simulate3b

apply(rbind(beta31[c(1:3,7,8,15,18,19),],beta32[c(4,5,6,10:13,16,20),]),2,mean) ## simulate3a
apply(rbind(beta32[c(1:3,7,8,15,18,19),],beta31[c(4,5,6,10:13,16,20),]),2,mean) ## simulate3b


c(alpha11[c(1:3,7,8,15,18,19)],alpha12[c(4,5,6,10:13,16,20)])
c(alpha12[c(1:3,7,8,15,18,19)],alpha11[c(4,5,6,10:13,16,20)])
c(alpha21[c(1:3,7,8,15,18,19)],alpha22[c(4,5,6,10:13,16,20)])
c(alpha22[c(1:3,7,8,15,18,19)],alpha21[c(4,5,6,10:13,16,20)])



bplot.sim <- function(beta_par, bTrue, name_title) {
  boxplot(beta_par, main = name_title)
  # create data for segments
  n <- ncol(beta_par)
  # width of each boxplot is 0.8
  x0s <- 1:n - 0.4
  x1s <- 1:n + 0.4
  # these are the y-coordinates for the horizontal lines
  # that you need to set to the desired values.
  y0s <- bTrue
  # add segments
  segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "red", lwd = 2)
}

bplot.sim(rbind(beta11[c(1:3,7,8,15,18,19),], beta12[c(4,5,6,10:13,16,20),]), 
          bTrue = c(-1.04674145, -0.70988616, 0.82395339), name = "beta11")
bplot.sim(rbind(beta12[c(1:3,7,8,15,18,19),], beta11[c(4,5,6,10:13,16,20),]), 
          bTrue = c(1.04674145, 0.70988616, -0.82395339), name = "beta12")

bplot.sim(rbind(beta21[c(1:3,7,8,15,18,19),],beta22[c(4,5,6,10:13,16,20),]), 
          bTrue = c(0.36239039, 0.50113490, 0.49765166), name = "beta21")
bplot.sim(rbind(beta22[c(1:3,7,8,15,18,19),],beta21[c(4,5,6,10:13,16,20),]), 
          bTrue = c(-1.36239039, 0.50113490, -0.29765166), name = "beta22")

bplot.sim(rbind(beta31[c(1:3,7,8,15,18,19),],beta22[c(4,5,6,10:13,16,20),]), 
          bTrue = c(-5.2947848, -0.9861937, 0.9005431), name = "beta31")
bplot.sim(rbind(beta32[c(1:3,7,8,15,18,19),],beta21[c(4,5,6,10:13,16,20),]), 
          bTrue = c(5.2947848, 0.9861937, -0.9005431), name = "beta32")
