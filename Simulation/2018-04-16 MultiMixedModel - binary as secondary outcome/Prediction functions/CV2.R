rm(list=ls(all=TRUE))

setwd("G:/2017-01 ROI/2017-01-03 Andrinopoulou - Uveitis project/Simulation/2018-04-16 MiltiMixedModel 2or3 out binary/Prediction functions")

##### packages
require(JMbayes)
require(rjags)
require(nlme)
require(lme4)
require(splines)


#################################
############## DATA #############
#################################
source("simulate2.R")



###### Sources
source(file="prediction functions2CV.R")


###### function for cross validation for the model
CVmodel <- function(ind, ttimes) {
  
  ###### training (for fitting) / testing (for prediction)
  testingdata <- data[data$id %in% ind & data$year <= ttimes, ]
  
  trainingdata <- data
  trainingdata <- trainingdata[trainingdata$id %in% ind & trainingdata$year <= ttimes, ]
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  test <- data[data$id %!in% ind, ]
  
  trainingdata <- rbind(test, trainingdata)  
  
  
  trainingdata <- trainingdata[complete.cases(trainingdata$year, trainingdata$ y3, trainingdata$drug, 
                                              trainingdata$y1, trainingdata$ y1), ]
  
  
  
  ###### fit the model
  source("analysis2.R", local = T)
  
  model.fit <- jags.model(file = "multiMixed2", data = Data, n.chains = con$n.chains, 
                          n.adapt = con$n.adapt, quiet = con$quiet)
  
  update(model.fit, con$n.burnin)
  res <- coda.samples(model.fit, parms,  n.iter = con$n.iter - con$n.burnin, thin = con$n.thin)
  codaFit <- as.mcmc.list(res)

  
  bss <- do.call(rbind,codaFit)
  n.sims <- nrow(bss)
  sims.list <- vector("list", length(parms))
  names(sims.list) <- parms
  for (p in seq_along(parms)) {
    ii <- grep(paste("^", parms[p], sep = ""), colnames(bss))
    sims.list[[p]] <- bss[, ii]
  }
  
  
  
  ###### Get predictions
  patdata <- list()
  ind <- unique(testingdata$id)
  for (k in 1:length(ind)){
    patdata[[k]] <- testingdata[testingdata$id == ind[k], ]
  }
  
  pred <- list()
  for (k in 1:length(ind)){
    ftimes <- data$year[data$id == ind[k]]
    if(  any(ftimes > ttimes) & !any(is.na(patdata[[k]]$y3 == T)) & !any(is.na(patdata[[k]]$year) == T) & !any(is.na(patdata[[k]]$drug) == T)  ){
      testtimes <- max(testingdata$year[testingdata$id == ind[k]])
      stimes <- unique(ftimes[ftimes > testtimes])
      stimes <- if (length(stimes) > 1) {
        sort(stimes, partial = length(stimes) - 1)[1:length(stimes)]
      } else if (length(stimes) == 1) {
        (ftimes[ftimes > testtimes])
      }
      stimes <- stimes[!is.na(stimes)]
      
      
      (pred[[k]] <- predictY(anadata = data, sims.list = sims.list, fm_VA = fm_VA, fm_I = fm_I, newdata = patdata[[k]], 
                             timesPred = stimes, 
                             M = 200L, timeVar = "year", idVar = "ideye", seed = 1L))
    }                     
  }
  
  #try(pred <- lapply(patdata, function(x) getPredictions(x, fitModel, horizon = ftimes, nsim = nSim,
  #                               setseed = 1, imaxit = 5, addError = F)$predData))
  
  #Warning message: 1: In e$fun(obj, substitute(ex), parent.frame(), e$data) : already exporting variable(s): patdata
  #save(pred, file = "pred.RData")
  #load("pred.RData")
  
  res <- list()
  for (k in 1:length(pred)){
    if (!is.null(pred[[k]])) {
      predData <- as.data.frame(pred[[k]])
      ftimes <- data$year[data$id == ind[k]]
      testtimes <- max(testingdata$year[testingdata$id == ind[k]])
      stimes <- unique(ftimes[ftimes > testtimes])
      stimes <- if (length(stimes) > 1) {
        sort(stimes, partial = length(stimes) - 1)[1:length(stimes)]
      } else if (length(stimes) == 1) {
        (ftimes[ftimes > testtimes])
      }
      stimes <- stimes[!is.na(stimes)]
      
      sumPat <- data[data$id == ind[k], ]
      vaReal <- sumPat[sumPat$year %in% c(stimes), c("id", "year", "y1", "y2")]
      colnames(vaReal) <- c("ideye", "year", "y1", "y2")
      
      vaPred <- predData[predData$times.to.pred %in% c(stimes), c("ideye", "times.to.pred", "ypred1", "ypredIN")] 
      colnames(vaPred) <- c("ideye", "year", "yPred", "INpred")
      dat <- merge(vaPred, vaReal, by = c("ideye", "year"))
      dat$ttime <- testtimes
      #if (dat$vaPred < 0) dat$vaPred <- 0
      dat$diff <-  if (any(dat$vaPred == 0) | any(dat$y1 == 0)){
        1*(((dat$yPred + 0.001) - (dat$y1 + 0.001)))
      } else {
        1*(((dat$yPred) - (dat$y1)))
      }
      res[[k]] <- dat
    }
  }
  
  do.call(rbind, res)
  #CVres1predcom <- do.call(rbind, res)
}
###########################################################



M = 2
CVres1pred <- list()
CVres2pred <- list()
CVres3pred <- list()


for (k in 1:M) {
  
  print(k)
  set.seed(k)
  V <- 5
  n <- nrow(data[!duplicated(data$id), ])
  splits <- split(seq_len(n), sample(rep(seq_len(V), length.out = n)))
  
  CVres1pred[[k]] <- lapply(splits, function(x) CVmodel(x, ttimes = 1))
  CVres2pred[[k]] <- lapply(splits, function(x) CVmodel(x, ttimes = 2))
  CVres3pred[[k]] <- lapply(splits, function(x) CVmodel(x, ttimes = 3))

  #fileout <- paste0("Results_modelNEW/resultsCVmodel", k,".RData")
  #save(CVres1pred, CVres2pred, CVres3pred, file = fileout)
  
}


CVres1predF <- list()
CVres2predF <- list()
CVres3predF <- list()

for (k in 1:M){
  
  #fileout <- paste0("Results_modelNEW/resultsCVmodel", k,".RData")
  #load(fileout)
  CVres1predF[[k]] <- do.call(rbind, CVres1pred[[k]])
  CVres2predF[[k]] <- do.call(rbind, CVres2pred[[k]])
  CVres3predF[[k]] <- do.call(rbind, CVres3pred[[k]])
  
}

CVres1predcom <- do.call(rbind, CVres1predF)

OkPred <- function(data, range){
  sum(abs(data$diff) <= range, na.rm = T)/length(data$diff)
}

OkPredoverall101 <- OkPred(CVres1predcom, 0.1)
