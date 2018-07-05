
n = 350 #n <- dim(mat1)[1]

K <- 10  # number of planned repeated measurements per subject, per outcome
t.max <- 19.5 # maximum follow-up time

################################################
# parameters for the linear mixed effects model 1
betaa1 <- c("(Intercept)" = 0.04674145, "Group1" = -1.70988616, "Time1" =  -0.82395339)
sigma1.y <- 1.017671  # measurement error standard deviation

# parameters for the logistic mixed effects model 
betaa2 <- c("(Intercept)" = -0.3151882, "Group1" = -0.6390350 , "Time1" =  -0.4082272)

# parameters for the logistic mixed effects model 
betaa3 <- c("(Intercept)" = -2.2947848, "Group1" = 1.9861937, "Time1" =  0.2005431)

# association parameter
alphaa1 <- 2.353078#-0.5 # association parameter - value
alphaa2 <- -0.9085359#-0.5 # association parameter - value


#D <- diag(c(6.9704385, 0.504228129, 1.52211686, 1.063437272))

invDvec <- c( 0.9335740,  0.7810344 , 2.0630577 , 0.8258380, -0.8279064, -0.3603033 , 
              0.7810344,  2.9226229 , 1.1579130 , 0.2087516, -0.8010016, -0.3315249 ,
              2.0630577,  1.1579130 , 6.2353232 , 2.8939485, -1.6217820 ,-0.7249335 ,
              0.8258380 , 0.2087516 , 2.8939485 , 4.8074160, -0.6370743 ,-1.2001506 ,
              -0.8279064, -0.8010016, -1.6217820 ,-0.6370743,  0.8903243 , 0.3994297 ,
              -0.3603033, -0.3315249 ,-0.7249335, -1.2001506 , 0.3994297 , 0.4546367 )

invDmat <- matrix(invDvec, 6, 6)
D <- solve(invDmat)

################################################

# design matrices for the longitudinal measurement model
# but this can be easily generalized
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max))))) # at which time points longitudinal measurements are supposed to be taken
group <- sample(c(0,1), n, replace = TRUE)

DF <- data.frame(year = times, drug = factor(rep(group, each = K)))

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

eta.y1 <- as.vector(X1 %*% betaa1 + rowSums(Z1 * b[id, 1:2])) + alphaa1 * y2 + alphaa2 * y3 # linear predictor
y1 <- rnorm(n * K, eta.y1, sigma1.y)


dat <- DF[, ]
dat$id <- id
dat$y1 <- y1
dat$y2 <- y2
dat$y3 <- y3

names(dat) <- c("year", "drug", "id", "y1", "y2", "y3")

dat.id <- dat[tapply(row.names(dat), dat$id, tail, 1), ]

data <- dat
data.id <- dat.id

