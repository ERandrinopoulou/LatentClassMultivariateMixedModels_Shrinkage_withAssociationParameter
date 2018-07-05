
#rm(list=ls(all=TRUE))

n = 350 #n <- dim(mat1)[1]
age <- rnorm(n, 45, 15.69578) #age <- mat1$age


##############################
#n <- 100 # number of subjects
##############################
K <- 10  # number of planned repeated measurements per subject, per outcome
t.max <- 19.5 # maximum follow-up time

################################################

# parameters for the linear mixed effects model 1
betas1 <- c("(Intercept)" = 8.0317, "Group1" = -5.8568, "Time1" =  -0.1578)
sigma1.y <- 0.6920  # measurement error standard deviation

# parameters for the linear mixed effects model 
betas2 <- c("(Intercept)" = -8.0317, "Group1" = 12.8568, "Time1" =  0.2)

# association parameter
alpha <-    0.3833#-0.5 # association parameter - value


D <- diag(c(0.9337, 0.1560, 0.8, 0.9)^2)

################################################

# design matrices for the longitudinal measurement model
# but this can be easily generalized
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max))))) # at which time points longitudinal measurements are supposed to be taken
group <- sample(c(0,1), n, replace = TRUE) #rep(0:1, each = n/2) # group indicator, i.e., '0' placebo, '1' active treatment
####################################
#age <- rnorm(n, 46.39673, 2.69578)
####################################
DF <- data.frame(year = times, drug = factor(rep(group, each = K)), age = rep(age, each = K))


#X <- model.matrix(~ 0 + drug  + ns(year, knots = kn, Boundary.knots = Bkn), data = DF)
#Z <- model.matrix(~ ns(year, knots = kn, Boundary.knots = Bkn), data = DF)
X1 <- model.matrix(~ drug + year, data = DF)
Z1 <- model.matrix(~ year, data = DF)

X2 <- model.matrix(~ drug + year, data = DF)
Z2 <- model.matrix(~ year, data = DF)

# design matrix for the survival model
W <- cbind("(Intercept)" = 1, "Age" = age)

################################################

#simulate random effects
library(MASS)

b <- mvrnorm(n, rep(0, nrow(D)), D)


# simulate longitudinal responses
id <- rep(1:n, each = K)

eta.y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b[id, 3:4])) # linear predictor
y2 <- rbinom(n * K, size = 1, prob = plogis(eta.y2))

eta.y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b[id, 1:2])) # linear predictor
y1 <- rnorm(n * K, eta.y1, sigma1.y)



dat <- DF[, ]
dat$id <- id
dat$y1 <- y1
dat$y2 <- y2

names(dat) <- c("year", "drug", "age", "id", "y1", "y2")

dat.id <- dat[tapply(row.names(dat), dat$id, tail, 1), ]


data <- dat
data.id <- dat.id

