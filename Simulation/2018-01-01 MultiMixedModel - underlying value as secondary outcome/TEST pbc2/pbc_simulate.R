
n = 350 #n <- dim(mat1)[1]

K <- 10  # number of planned repeated measurements per subject, per outcome
t.max <- 19.5 # maximum follow-up time

################################################
# parameters for the linear mixed effects model 1
betaa1 <- c("(Intercept)" = 6.2835566, "Group1" = -0.1019855, "Time1" =  -0.4281292)
sigma1.y <- 2.815901  # measurement error standard deviation

# parameters for the linear mixed effects model 
betaa2 <- c("(Intercept)" = 3.1006302, "Group1" = -0.6101297, "Time1" =  -0.3712205)

# association parameter
alphaa <- 0.2820991#-0.5 # association parameter - value


#D <- diag(c(6.9704385, 0.504228129, 1.52211686, 1.063437272))

invDvec <- c (0.29876080 , 2.61208478, -0.02967438 ,-0.06063465,
             2.61208478 ,24.78371981 ,-0.12133228 ,-1.40824154,
             -0.02967438, -0.12133228 , 0.33281734, 1.45271670,
             -0.06063465, -1.40824154,  1.45271670 ,21.41737861 )

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

eta.y1 <- as.vector(X1 %*% betaa1 + rowSums(Z1 * b[id, 1:2])) + alphaa * eta.y2  # linear predictor
y1 <- rnorm(n * K, eta.y1, sigma1.y)


pbc2 <- DF[, ]
pbc2$id <- id
pbc2$y1 <- y1
pbc2$y2 <- y2

names(pbc2) <- c("years", "drug", "id", "serBilir", "hepatomegaly")



