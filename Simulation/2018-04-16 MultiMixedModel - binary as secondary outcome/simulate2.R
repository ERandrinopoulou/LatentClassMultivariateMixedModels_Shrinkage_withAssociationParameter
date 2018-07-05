
n = 350 #n <- dim(mat1)[1]

K <- 10  # number of planned repeated measurements per subject, per outcome
t.max <- 19.5 # maximum follow-up time

################################################
# parameters for the linear mixed effects model 1
betaa1 <- c("(Intercept)" = 0.04674145, "Group1" = -1.70988616, "Time1" =  -0.82395339)
sigma1.y <- 1.017671  # measurement error standard deviation

# parameters for the linear mixed effects model 
betaa2 <- c("(Intercept)" = -0.36239039, "Group1" = -0.50113490, "Time1" =  -0.39765166)

# association parameter
alphaa <- -2.074567#-0.5 # association parameter - value


#D <- diag(c(6.9704385, 0.504228129, 1.52211686, 1.063437272))

invDvec <- c(0.9455738,  0.7497205, -2.135328, -0.9877032,
             0.7497205,  2.6045574, -1.847133, -0.8914816,
             -2.1353280, -1.8471326,  5.713290,  2.7678937,
             -0.9877032, -0.8914816,  2.767894,  2.2962598)

invDmat <- matrix(invDvec, 4, 4)
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

################################################

#simulate random effects
library(MASS)

b <- mvrnorm(n, rep(0, nrow(D)), D)


# simulate longitudinal responses
id <- rep(1:n, each = K)

eta.y2 <- as.vector(X2 %*% betaa2 + rowSums(Z2 * b[id, 3:4])) # linear predictor
y2 <- rbinom(n * K, size = 1, prob = plogis(eta.y2))

eta.y1 <- as.vector(X1 %*% betaa1 + rowSums(Z1 * b[id, 1:2])) + alphaa * y2  # linear predictor
y1 <- rnorm(n * K, eta.y1, sigma1.y)


dat <- DF[, ]
dat$id <- id
dat$y1 <- y1
dat$y2 <- y2

names(dat) <- c("year", "drug", "id", "y1", "y2")

dat.id <- dat[tapply(row.names(dat), dat$id, tail, 1), ]

data <- dat
data.id <- dat.id

