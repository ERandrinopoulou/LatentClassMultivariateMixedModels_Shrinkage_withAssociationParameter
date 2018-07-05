

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



mu01 <- rep(0,(ncZ+ncZ2+ncZ3))
mu02 <- rep(0,(ncZ+ncZ2+ncZ3))




betas11 <- rep(0, ncX)
var.betas11 <- rep(con$K, ncX)
betas12 <- rep(0, ncX)
var.betas12 <- rep(con$K, ncX)


betas21 <- rep(0, ncX2)
var.betas21 <- rep(con$K, ncX2)
betas22 <- rep(0, ncX2)
var.betas22 <- rep(con$K, ncX2)


betas31 <- rep(0, ncX2)
var.betas31 <- rep(con$K, ncX2)
betas32 <- rep(0, ncX2)
var.betas32 <- rep(con$K, ncX2)



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

             
             mu01 = mu01, 
             mu02 = mu02, 
             
             prior.cl = rep(15.5, 2),
             
             priorMean.betas11 = betas11, 
             priorTau.betas11 = diag(1/var.betas11),
             priorMean.betas12 = betas12, 
             priorTau.betas12 = diag(1/var.betas12),
             
             priorMean.betas21 = betas21, 
             priorTau.betas21 = diag(1/var.betas21),
             priorMean.betas22 = betas22, 
             priorTau.betas22 = diag(1/var.betas22),
             
             priorMean.betas31 = betas31, 
             priorTau.betas31 = diag(1/var.betas31),
             priorMean.betas32 = betas32, 
             priorTau.betas32 = diag(1/var.betas32),
             
             priorA.tau = (1/sigma2)^2/10,
             priorB.tau = (1/sigma2)/10, 
             
             priorR.D1 = diag(1,(ncZ+ncZ2+ncZ3)), priorK.D1 = (ncZ+ncZ2+ncZ3),
             priorR.D2 = diag(1,(ncZ+ncZ+ncZ3)), priorK.D2 = (ncZ+ncZ2+ncZ3)
             
)

parms <- c("betas11", "betas12", "betas21", "betas22", "betas31", "betas32", "alphas11", "alphas12", "alphas21", "alphas22", "tau", "inv.D1", "inv.D2", "pr")
#####################



model.fit <- jags.model(file = "multiMixed3shrink", data = Data, n.chains = con$n.chains, 
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




betas11 <- apply(sims.list[[1]][,1:(ncX)], 2, mean)  # betas1 
betas12 <- apply(sims.list[[2]][,1:(ncX)], 2, mean)  # betas1 


betas21 <- apply(sims.list[[3]][,1:ncX2], 2, mean)  # betas2 
betas22 <- apply(sims.list[[4]][,1:ncX2], 2, mean)  # betas2 

betas31 <- apply(sims.list[[5]][,1:ncX2], 2, mean)  # betas2 
betas32 <- apply(sims.list[[6]][,1:ncX2], 2, mean)  # betas2 


alphas11 <- mean(sims.list[[7]])          # alpha11
alphas12 <- mean(sims.list[[8]])          # alpha12
alphas21 <- mean(sims.list[[9]])          # alpha21
alphas22 <- mean(sims.list[[10]])         # alpha22


tau <- mean(sims.list[[11]])          # tau


invD1 <- apply(sims.list[[12]],2,mean)  # invD
invD2 <- apply(sims.list[[13]],2,mean)  # invD

prr <- sims.list$pr

