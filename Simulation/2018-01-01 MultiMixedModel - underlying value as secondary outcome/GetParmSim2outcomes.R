
rm(list=ls(all=TRUE))

setwd("C:/Users/015004/Desktop/2017-01 ROI/2017-01-03/JM_anaPred underlying value")

##### packages
require(JMbayes)
require(rjags)
require(nlme)
require(lme4)
require(splines)


###### Sources


anadata <-  read.table("C:/Users/015004/Desktop/2017-01 ROI/2017-01-03/01anadata.csv", header = TRUE, sep = ",", row.names= 1)


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


############################


fm_VA <- lme(y ~ ageelder + time, data = anadata, na.action = na.exclude, 
             random =~ time|ideye)


lmeObject <- fm_VA

anadata$anysurgery <- as.factor(anadata$anysurgery)
anadata$sc <- as.factor(anadata$sc)

fm_I <- glmer(x ~ ageelder + time + (time | ideye), data = anadata, 
              family = binomial, 
              na.action = na.exclude,  control=glmerControl(optimizer = "bobyqa",
                                                            optCtrl = list(maxfun = 100000)))


lmeObject2 <- fm_I


timeVar <- "time"
lag <- 0

#####################################################

######################################################


id <- as.numeric(anadata$ideye) #as.vector(unclass(lmeObject$groups[[1]]))


offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))

#################
formYx <- formula(fm_VA)
TermsX <- lmeObject$terms
mfX <- model.frame(TermsX, data = anadata)
X <- model.matrix(formYx, mfX)

formYz <- formula(fm_VA$modelStruct$reStruct[[1]])
mfZ <- model.frame(terms(formYz), data = anadata)
TermsZ <- attr(mfZ, "terms")
Z <- model.matrix(formYz, mfZ)

#


###
formYx2 <- formula(noquote(paste("x ~", c(terms(fm_I)[[3]])))) #formula(lmeObject2)
TermsX2 <-  terms(lmeObject2)
mfX2 <- model.frame(TermsX2, data = anadata)
X2 <- model.matrix(formYx2, mfX2)

formYz2 <- formula(noquote(paste("~", attributes(terms(fm_I))$term.labels[[1]]))) #~ ns(time, 3) #time
mfZ2 <- model.frame(terms(formYz2), data = anadata)
TermsZ2 <- attr(mfZ2, "terms")
Z2 <- model.matrix(formYz2, mfZ2)

####################################################

y.long <- model.response(mfX, "numeric")
y.long2 <- model.response(mfX2, "numeric")
y <- list(y = y.long, offset = offset, y2 = y.long2, 
          lag = lag)



####################



con <- list(program = "JAGS", n.chains = 1, n.iter = 250000,
            n.burnin = 125000, n.thin = 10, n.adapt = 1500, K = 100,
            C = 5000, working.directory = getwd(), bugs.directory = "C:/Program Files/WinBUGS14/",
            openbugs.directory = NULL, clearWD = TRUE, over.relax = TRUE,
            knots = NULL, ObsTimes.knots = TRUE, lng.in.kn = 5, ordSpline = 4,
            bugs.seed = 1, quiet = FALSE)




#################################
ncX <- ncol(X)
ncZ <- ncol(Z)

ncX2 <- ncol(X2)
ncZ2 <- ncol(Z2)



nb <- ncZ + ncZ2


mu0 <- rep(0,(ncZ+ncZ2))


betas <- rep(0, ncX)
var.betas <- rep(con$K, ncX)

betas2 <- rep(0, ncX2)
var.betas2 <- rep(con$K, ncX2)

alphas <- 0
var.alphas <- con$K

b <- cbind(data.matrix(random.effects(lmeObject)),data.matrix(ranef(lmeObject2)$id))


nY <- nrow(b)
sigma2 <- lmeObject$sigma^2



#####################################################################################################

Data <- list(N = nY, offset = offset, X = X, 
             y = y$y, 
             Z = Z, 
             
             
             X2 = X2, 
             y2 = y$y2, 
             Z2 = Z2, 
             nb = nb, 
             ncX = ncol(X), 
             ncZ = ncol(Z), 
             ncX2 = ncol(X2),
             ncZ2 = ncol(Z2), 
             
             
             mu0 = mu0, 
             priorMean.betas = betas, 
             priorTau.betas = diag(1/var.betas),
             
             priorMean.betas2 = betas2, 
             priorTau.betas2 = diag(1/var.betas2),
             
             priorA.tau = (1/sigma2)^2/10,
             priorB.tau = (1/sigma2)/10, 
             
             priorMean.alphas = alphas,
             priorTau.alphas = 1/var.alphas,
             
             priorR.D = diag(1,(ncZ+ncZ2)), priorK.D = (ncZ+ncZ2)
             
)

parms <- c("betas","betas2", "tau", "inv.D", "b", "alphas")
#####################
model.fit <- jags.model(file = "multiMixed.txt", data = Data, n.chains = con$n.chains, 
                        n.adapt = con$n.adapt, quiet = con$quiet)

update(model.fit, con$n.burnin)
res <- coda.samples(model.fit, parms,  n.iter = con$n.iter - con$n.burnin, thin = con$n.thin)
codaFit <- as.mcmc.list(res)

######################################################################################
###################################### OUTPUT ########################################
######################################################################################

bss <- do.call(rbind,codaFit)
n.sims <- nrow(bss)
sims.list <- vector("list", length(parms))
names(sims.list) <- parms
for (p in seq_along(parms)) {
  ii <- grep(paste("^", parms[p], sep = ""), colnames(bss))
  sims.list[[p]] <- bss[, ii]
}

stdErr <-
  function (x) {
    x <- as.matrix(x)
    vars <- apply(x, 2L, var)
    ess <- effectiveSize(x)
    sqrt(vars / ess)
  }

computeP <-
  function (x) {
    above <- mean(x >= 0)
    below <- mean(x < 0)
    2 * min(above, below)
  }




indb <- (parms) != "b"
PostMean = lapply(sims.list[indb], function (x) apply(as.matrix(x), 2L, mean))
PostSE = lapply(sims.list[indb], function (x) apply(as.matrix(x), 2L, stdErr))
Pvalues = lapply(sims.list[indb], function (x) apply(as.matrix(x), 2L, computeP))



################ effect plot #######################

newdata <- with(anadata, data.frame(
  time = rep(seq(min(anadata$time, max(anadata$time)), length = 70), times = 2),
  complication = rep(0, 140),
  age = rep(median(anadata$age), 140),
  IN = rep(c(0,1), each = 70)
))


newdata$complication <- factor(newdata$complication, levels = c(0, 1), labels = c("yes", "no")) 
newdata$IN <- factor(newdata$IN, levels = c(0, 1), labels = c("IN yes", "IN no")) 

X <- model.matrix(~ ns(time, 5) + complication + age + IN, data = newdata)


newdata$pred <- X %*% c(PostMean$betas[1:8], PostMean$alphas)
loCI <- apply(cbind(sims.list$betas[,1:8], sims.list$alphas),2, function(x) { quantile(x, c(0.025))} )
upCI  <- apply(cbind(sims.list$betas[,1:8], sims.list$alphas),2, function(x) { quantile(x, c(0.975))} )
newdata$lo <-  X %*% loCI
newdata$up <-  X %*% upCI

backtransva <- function( y , epsilon1 = 1e-03 , epsilon2 = 1e-01 )
{
  ( ( 1 + epsilon2 ) * exp( y ) - epsilon1 )/ ( 1 + exp( y ) ) 
}



xyplot(backtransva(pred) + backtransva(lo) + backtransva(up) ~ time | IN, data = newdata, type = "l", lty = c(1,2,2),
       col = c(1,1,1), lwd = 2, ylab = "Visual acuity", xlab = "Time (years)")

