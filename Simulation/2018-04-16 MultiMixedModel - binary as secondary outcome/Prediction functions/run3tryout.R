rm(list=ls(all=TRUE))


setwd("G:/2017-01 ROI/2017-01-03 Andrinopoulou - Uveitis project/Simulation/2018-04-16 MiltiMixedModel 2or3 out binary/Prediction functions")

##### packages
require(JMbayes)
require(rjags)
require(nlme)
require(lme4)
require(splines)
require(lattice)

##################################################
##################### DATA #######################
##################################################
anadata <-  read.table("G:/2017-01 ROI/2017-01-03 Andrinopoulou - Uveitis project/Simulation/2018-04-16 MiltiMixedModel 2or3 out binary/Prediction functions/01anadata.csv", header = TRUE, sep = ",", row.names= 1)

###### Prepare data
anadata$Age <- anadata$age

anadata <- within(anadata,{
  age <- (Age - 43)/10
  age2 <- age^2 
  ageyoung <- 1*(Age<=18)
  ageelder <- 1*(Age>=60)
  sc <- ifelse(startcom<=6*7,0,1)
  treat <- factor(treat)
  spellf <- factor(spell)
})

anadata <- anadata[complete.cases(anadata$y, anadata$time, anadata$complication, anadata$age, anadata$x, anadata$treat, anadata$anysurgery,
                                  anadata$sc), ]

##################################################
################### ANALYSIS #####################
##################################################
#source("Analysis.R")


load("resSim3.RData")



##################################################
############## EXPORT RESULTS ####################
##################################################
parms <- c("betas", "betas2", "betas3", "tau", "inv.D", "b")

bss <- do.call(rbind,codaFit)

n.sims <- nrow(bss)
sims.list <- vector("list", length(parms))
names(sims.list) <- parms
for (p in seq_along(parms)) {
  ii <- grep(paste("^", parms[p], sep = ""), colnames(bss))
  sims.list[[p]] <- bss[, ii]
}

##################################################
############# RUN PREDICTIONS ####################
##################################################
source(file="prediction functions3.R")


##################################################
################ RUN MODELS ######################
##################################################

anadata$treat[anadata$treat == 2] <- 1
anadata$treat <- droplevels(anadata$treat)

fm_VA <- lme(y ~ treat + time, data = anadata, na.action = na.exclude, 
             random =~ time|ideye, control=list(maxIter=1000, 
                                              msMaxIter=1000, niterEM=1000))



fm_I <- glmer(x ~ treat + time + (time | ideye), data = anadata, 
              family = binomial, 
              na.action = na.exclude,  control=glmerControl(optimizer = "bobyqa",
                                                            optCtrl = list(maxfun = 100000)))


fm_Com <- glmer(complication ~ treat + time + (time | ideye), data = anadata, 
                family = binomial, 
                na.action = na.exclude,  control=glmerControl(optimizer = "bobyqa",
                                                              optCtrl = list(maxfun = 100000)))

##################################################
############## RUN PREDICTIONS ###################
##################################################

indpred <- "116OD"
indpred <- "116OS"
indpred <- "268OD"
indpred <- "114OS"
indpred <- "136OD"
indpred <- "32OD"
newdata <- anadata[anadata$ideye %in% indpred, ]

i=3
stimes <- newdata$time[newdata$time > newdata$time[i]]
stimes <- seq(min(stimes), max(stimes), length.out = 20)
newdata1 <- newdata[newdata$time <= min(stimes),]

anadata = anadata
sims.list = sims.list
fm_VA = fm_VA
fm_I = fm_I
newdata = newdata1
timesPred = stimes
M = 200L
timeVar = "time"
idVar = "ideye"
seed = 1L


(pred <- predictY(anadata = anadata, sims.list = sims.list, fm_VA = fm_VA, fm_I = fm_I, newdata = patdata[[k]], timesPred = stimes, 
                       M = 200L, timeVar = "time", idVar = "ideye", seed = 1L))
