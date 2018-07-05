
#rm(list=ls(all=TRUE))

n = 1050#350 #n <- dim(mat1)[1]
age <- rnorm(n, 45, 15.69578) #age <- mat1$age


##############################
#n <- 100 # number of subjects
##############################
K <- 10  # number of planned repeated measurements per subject, per outcome
t.max <- 19.5 # maximum follow-up time

################################################

# parameters for the linear mixed effects model 1
betaa1 <- c("(Intercept)" = 5.0317, "Group1" = -5.8568, "Time1" =  -0.1578)
sigma1.y <- 0.6920  # measurement error standard deviation

# parameters for the logistic mixed effects model 
betaa2 <- c("(Intercept)" = -5.0317, "Group1" = 5.8568, "Time1" =  0.2)

# parameters for the logistic mixed effects model 
betaa3 <- c("(Intercept)" = 1.0317, "Group1" = 0.8568, "Time1" =  0.5)


# association parameters
alphaa1 <- 0.3833 # association parameter - value
alphaa2 <- 0.03 # association parameter - value


D <- diag(c(0.9337, 0.1560, 0.9337, 0.5560, 0.9337, 0.5560)^2)

################################################

# design matrices for the longitudinal measurement model
# but this can be easily generalized
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max))))) # at which time points longitudinal measurements are supposed to be taken
group <- sample(c(0,1), n, replace = TRUE) #rep(0:1, each = n/2) # group indicator, i.e., '0' placebo, '1' active treatment
####################################
#age <- rnorm(n, 46.39673, 2.69578)
####################################
DF <- data.frame(year = times, drug = factor(rep(group, each = K)), age = rep(age, each = K))


X1 <- model.matrix(~ drug + year, data = DF)
Z1 <- model.matrix(~ year, data = DF)

X2 <- model.matrix(~ drug + year, data = DF)
Z2 <- model.matrix(~ year, data = DF)

X3 <- model.matrix(~ drug + year, data = DF)
Z3 <- model.matrix(~ year, data = DF)

################################################

#simulate random effects
library(MASS)

b <- mvrnorm(n, rep(0, nrow(D)), D)


# simulate longitudinal responses
id <- rep(1:n, each = K)


eta.y2 <- as.vector(X2 %*% betaa2 + rowSums(Z2 * b[id, 3:4])) # linear predictor
y2 <- rbinom(n * K, size = 1, prob = plogis(eta.y2))

eta.y3 <- as.vector(X3 %*% betaa3 + rowSums(Z3 * b[id, 5:6])) # linear predictor
y3 <- rbinom(n * K, size = 1, prob = plogis(eta.y3))

eta.y1 <- as.vector(X1 %*% betaa1 + rowSums(Z1 * b[id, 1:2]) + alphaa1*eta.y2 + alphaa2*eta.y3) # linear predictor
y1 <- rnorm(n * K, eta.y1, sigma1.y)



dat <- DF[, ]
dat$id <- id
dat$y1 <- y1
dat$y2 <- y2
dat$y3 <- y3

names(dat) <- c("year", "drug", "age", "id", "y1", "y2", "y3")

dat.id <- dat[tapply(row.names(dat), dat$id, tail, 1), ]


data <- dat
data.id <- dat.id

