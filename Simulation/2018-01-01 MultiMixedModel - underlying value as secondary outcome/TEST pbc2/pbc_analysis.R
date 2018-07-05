
pbc2 <- pbc2[order(pbc2$id, pbc2$years),]

#pbc2$hepatomegaly <- as.numeric(pbc2$hepatomegaly) - 1


pbc2 <- pbc2[complete.cases(pbc2$serBilir, pbc2$drug, pbc2$years, pbc2$hepatomegaly), ]

############################

lmeObject1 <- lme(serBilir ~ drug + years, data = pbc2, na.action = na.exclude, 
             random =~ years | id)


lmeObject2 <- glmer(hepatomegaly ~ drug + years + (years | id), data = pbc2, 
              family = binomial, control = glmerControl(optimizer = "bobyqa"),
              na.action = na.exclude)



timeVar <- "years"
lag <- 0

######################################################


id <- as.numeric(pbc2$id) #as.vector(unclass(lmeObject$groups[[1]]))


offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))

#################
formYx <- formula(lmeObject1)
TermsX <- lmeObject1$terms
mfX <- model.frame(TermsX, data = pbc2)
X <- model.matrix(formYx, mfX)

formYz <- formula(lmeObject1$modelStruct$reStruct[[1]])
mfZ <- model.frame(terms(formYz), data = pbc2)
TermsZ <- attr(mfZ, "terms")
Z <- model.matrix(formYz, mfZ)

#


###
formYx2 <- ~ drug + years #formula(lmeObject2)
TermsX2 <-  terms(lmeObject2)
mfX2 <- model.frame(TermsX2, data = pbc2)
X2 <- model.matrix(formYx2, mfX2)

formYz2 <- ~ years
mfZ2 <- model.frame(terms(formYz2), data = pbc2)
TermsZ2 <- attr(mfZ2, "terms")
Z2 <- model.matrix(formYz2, mfZ2)
####################################################

y.long <- model.response(mfX, "numeric")
y.long2 <- model.response(mfX2, "numeric")
y <- list(y = y.long, offset = offset, y2 = y.long2, 
          lag = lag)




####################


con <- list(program = "JAGS", n.chains = 1, n.iter = 25000,
            n.burnin = 20000, n.thin = 2, n.adapt = 1500, K = 100,
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


betas <- rep(0, ncX + 1)
var.betas <- rep(con$K, ncX + 1)

betas2 <- rep(0, ncX2)
var.betas2 <- rep(con$K, ncX2)


b <- cbind(data.matrix(random.effects(lmeObject1)),data.matrix(ranef(lmeObject2)$id))


nY <- nrow(b)
sigma2 <- lmeObject1$sigma^2

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
             
             priorR.D = diag(1,(ncZ+ncZ2)), priorK.D = (ncZ+ncZ2)
             
)

parms <- c("betas","betas2", "tau", "inv.D", "b")
#####################


model.fit <- jags.model(file = "multiMixed", data = Data, n.chains = con$n.chains, 
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




betas1 <- apply(sims.list[[1]][,1:(ncX + 1)],2,mean)  # betas1 


betas2 <- apply(sims.list[[2]][,1:ncX2],2,mean)  # betas2 



tau <- mean(sims.list[[3]])          # tau



invD <- apply(sims.list[[4]],2,mean)  # invD
