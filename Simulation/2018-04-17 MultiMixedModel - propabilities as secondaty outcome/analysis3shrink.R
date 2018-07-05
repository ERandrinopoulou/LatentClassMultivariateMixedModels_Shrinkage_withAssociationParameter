

############################

fm_VA <- lme(y1 ~ drug + year, data = data, na.action = na.exclude, 
             random =~ year | id)

lmeObject <- fm_VA

fm_I <- glmer(y2 ~ drug + year + (year | id), data = data, 
              family = binomial, control = glmerControl(optimizer = "bobyqa"),
              na.action = na.exclude)


lmeObject2 <- fm_I

fm_Com <- glmer(y3 ~ drug + year + (year | id), data = data, 
              family = binomial, control = glmerControl(optimizer = "bobyqa"),
              na.action = na.exclude)


lmeObject3 <- fm_Com

timeVar <- "year"
lag <- 0

#####################################################

######################################################


id <- as.numeric(data$id) #as.vector(unclass(lmeObject$groups[[1]]))


offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))

#################
formYx <- formula(lmeObject)
TermsX <- lmeObject$terms
mfX <- model.frame(TermsX, data = data)
X <- model.matrix(formYx, mfX)

formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
mfZ <- model.frame(terms(formYz), data = data)
TermsZ <- attr(mfZ, "terms")
Z <- model.matrix(formYz, mfZ)

#


###
formYx2 <- ~ drug + year #formula(lmeObject2)
TermsX2 <-  terms(lmeObject2)
mfX2 <- model.frame(TermsX2, data = data)
X2 <- model.matrix(formYx2, mfX2)

formYz2 <- ~ year
mfZ2 <- model.frame(terms(formYz2), data = data)
TermsZ2 <- attr(mfZ2, "terms")
Z2 <- model.matrix(formYz2, mfZ2)

###
formYx3 <- ~ drug + year #formula(lmeObject2)
TermsX3 <-  terms(lmeObject3)
mfX3 <- model.frame(TermsX3, data = data)
X3 <- model.matrix(formYx3, mfX3)

formYz3 <- ~ year
mfZ3 <- model.frame(terms(formYz3), data = data)
TermsZ3 <- attr(mfZ3, "terms")
Z3 <- model.matrix(formYz3, mfZ2)
####################################################

y.long <- model.response(mfX, "numeric")
y.long2 <- model.response(mfX2, "numeric")
y.long3 <- model.response(mfX3, "numeric")
y <- list(y = y.long, offset = offset, y2 = y.long2, y3 = y.long3,
          lag = lag)



####################



con <- list(program = "JAGS", n.chains = 1, n.iter = 150000,
            n.burnin = 100000, n.thin = 2, n.adapt = 1500, K = 100,
            C = 5000, working.directory = getwd(), bugs.directory = "C:/Program Files/WinBUGS14/",
            openbugs.directory = NULL, clearWD = TRUE, over.relax = TRUE,
            knots = NULL, ObsTimes.knots = TRUE, lng.in.kn = 5, ordSpline = 4,
            bugs.seed = 1, quiet = FALSE)


#################################
ncX <- ncol(X)
ncZ <- ncol(Z)

ncX2 <- ncol(X2)
ncZ2 <- ncol(Z2)

ncX3 <- ncol(X3)
ncZ3 <- ncol(Z3)

nb <- ncZ + ncZ2 + ncZ3


mu0 <- rep(0,(ncZ+ncZ2+ncZ3))


betas1 <- rep(0, ncX)
var.betas1 <- rep(con$K, ncX)

betas2 <- rep(0, ncX2)
var.betas2 <- rep(con$K, ncX2)


betas3 <- rep(0, ncX2)
var.betas3 <- rep(con$K, ncX2)



b <- cbind(data.matrix(random.effects(lmeObject)), data.matrix(ranef(lmeObject2)$id), data.matrix(ranef(lmeObject3)$id))


nY <- nrow(b)
sigma2 <- lmeObject$sigma^2

#####################################################################################################

Data <- list(N = nY, offset = offset, X = X, 
             y = y$y, 
             Z = Z, 
             
             
             X2 = X2, 
             y2 = y$y2, 
             Z2 = Z2, 
             X3 = X3, 
             y3 = y$y3, 
             Z3 = Z3, 
             
             nb = nb, 
             ncX = ncol(X), 
             ncZ = ncol(Z), 
             ncX2 = ncol(X2),
             ncZ2 = ncol(Z2), 
             ncX3 = ncol(X3),
             ncZ3 = ncol(Z3), 
             
             
             mu0 = mu0, 

             priorMean.betas1 = betas1, 
             priorTau.betas1 = diag(1/var.betas1),
             
             priorMean.betas2 = betas2, 
             priorTau.betas2 = diag(1/var.betas2),
             
             priorMean.betas3 = betas3, 
             priorTau.betas3 = diag(1/var.betas3),
    
             priorA.tau = (1/sigma2)^2/10,
             priorB.tau = (1/sigma2)/10, 
             
             priorR.D = diag(1,(ncZ+ncZ2+ncZ3)), priorK.D = (ncZ+ncZ2+ncZ3)
             
)

parms <- c("betas", "betas2", "betas3", alphas1", "alphas2","tau", "inv.D", "b")
#####################



model.fit <- jags.model(file = "multiMixed3shrink.R", data = Data, n.chains = con$n.chains, 
                        n.adapt = con$n.adapt, quiet = con$quiet)


update(model.fit, con$n.burnin)
res <- coda.samples(model.fit, parms,  n.iter = con$n.iter - con$n.burnin, thin = con$n.thin)
codaFit <- as.mcmc.list(res)



bss <- do.call(rbind,codaFit)
colnames(bss)
n.sims <- nrow(bss)
sims.list <- vector("list", length(parms))
names(sims.list) <- parms
for (p in seq_along(parms)) {
  ii <- grep(paste("^", parms[p], sep = ""), colnames(bss))
  sims.list[[p]] <- bss[, ii]
}




betas1 <- apply(sims.list[[1]][,1:(ncX + 2)],2,mean)  # betas1 


betas2 <- apply(sims.list[[2]][,1:ncX2],2,mean)  # betas2 

betas3 <- apply(sims.list[[3]][,1:ncX3],2,mean)  # betas2 


tau <- mean(sims.list[[4]])          # tau



invD <- apply(sims.list[[5]],2,mean)  # invD

