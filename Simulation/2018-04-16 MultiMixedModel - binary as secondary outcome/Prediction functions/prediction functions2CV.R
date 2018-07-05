# rm(list=ls(all=TRUE))
# 
# ##############################################
# ########## LIBRARIES #########################
# ##############################################
library("mvtnorm")
# 
# 
# ##################################################
# ############## LOAD MODEL ########################
# ##################################################
# load("F:/2017-01 ROI/2017-01-03/resMVGstanDimitris.RData")
# 

####################################################
################## FUNCTIONS #######################
####################################################
pii <- function(X, betas, Z, b) {
  exp(X%*%betas + Z%*%b)/(1 + exp(X%*%betas + Z%*%b)) #1/(1 + exp(-X%*%betas - Z%*%b))
}

like <- function(p, y){
  dbinom(y, 1, p, log = TRUE)
}


backtransva <- function( y , epsilon1 = 1e-03 , epsilon2 = 1e-01 )
{
  ( ( 1 + epsilon2 ) * exp( y ) - epsilon1 )/ ( 1 + exp( y ) ) 
}


cd <-
  function (x, f, ..., eps = .Machine$double.eps^0.25) {
    n <- length(x)
    res <- numeric(n)
    ex <- eps * (abs(x) + eps)
    for (i in seq_len(n)) {
      x1 <- x2 <- x
      x1[i] <- x[i] + ex[i]
      x2[i] <- x[i] - ex[i]
      diff.f <- c(f(x1, ...) - f(x2, ...))
      diff.x <- x1[i] - x2[i]
      res[i] <- diff.f / diff.x
    }
    res
  }

rmvt <-
  function (n, mu, Sigma, df) {
    p <- length(mu)
    if (is.list(Sigma)) {
      ev <- Sigma$values
      evec <- Sigma$vectors
    } else {
      ed <- eigen(Sigma, symmetric = TRUE)
      ev <- ed$values
      evec <- ed$vectors
    }
    X <- drop(mu) + tcrossprod(evec * rep(sqrt(pmax(ev, 0)), each = p), 
                               matrix(rnorm(n * p), n)) / rep(sqrt(rchisq(n, df)/df), each = p)
    if (n == 1L) drop(X) else t.default(X)
  }


dmvt <-
  function (x, mu, Sigma = NULL, invSigma = NULL, df, log = FALSE, prop = TRUE) {
    if (!is.numeric(x)) 
      stop("'x' must be a numeric matrix or vector")
    if (!is.matrix(x)) 
      x <- rbind(x)
    p <- length(mu)
    if (is.null(Sigma) && is.null(invSigma))
      stop("'Sigma' or 'invSigma' must be given.")
    if (!is.null(Sigma)) {
      if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
      } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors
      }
      if (!all(ev >= -1e-06 * abs(ev[1]))) 
        stop("'Sigma' is not positive definite")
      invSigma <- evec %*% (t(evec)/ev)
      if (!prop)
        logdetSigma <- sum(log(ev))
    } else {
      if (!prop)
        logdetSigma <- c(-determinant(invSigma)$modulus)
    }
    ss <- x - rep(mu, each = nrow(x))
    quad <- rowSums((ss %*% invSigma) * ss)/df
    if (!prop)
      fact <- lgamma((df + p)/2) - lgamma(df/2) - 0.5 * (p * (log(pi) + 
                                                                log(df)) + logdetSigma)
    if (log) {
      if (!prop) as.vector(fact - 0.5 * (df + p) * log(1 + quad)) else as.vector(- 0.5 * (df + p) * log(1 + quad))
    } else {
      if (!prop) as.vector(exp(fact) * ((1 + quad)^(-(df + p)/2))) else as.vector(((1 + quad)^(-(df + p)/2)))
    }
  }


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



predictY <- function(anadata, sims.list, fm_VA, fm_I, newdata, timesPred, M = 200L, 
                     timeVar = "time", idVar = "ideye", seed = 1L) {
  
  
  indpred <- newdata$id
  
  y1 <- newdata$y1[newdata$id %in% indpred]
  y2 <- newdata$y2[newdata$id %in% indpred]
  
  
  ########## Design matrices for newdata
  formYx_y <- formula(fm_VA)
  mfX_y <- model.frame(terms(formYx_y), data = anadata)
  TermsX1 <- attr(mfX_y, "terms")
  formYz_y <- formula(fm_VA$modelStruct$reStruct[[1]])#time #formula(lmeObject$modelStruct$reStruct[[1]])
  mfZ_y <- model.frame(terms(formYz_y), data = anadata)
  TermsZ1 <- attr(mfZ_y, "terms")
  formYx_x <- formula(noquote(paste("y2 ~", c(terms(fm_I)[[3]]))))
  mfX_x <- model.frame(terms(formYx_x), data = anadata)
  TermsX2 <- attr(mfX_x, "terms")
  formYz_x <- formula(noquote(paste("~", attributes(terms(fm_I))$term.labels[[1]]))) #~ ns(time, 3) #time
  mfZ_x <- model.frame(terms(formYz_x), data = anadata)
  TermsZ2 <- attr(mfZ_x, "terms")
  
  
  mfXpred1 <- model.frame(TermsX1, data = newdata)
  mfZpred1 <- model.frame(TermsZ1, data = newdata)
  mfXpred2 <- model.frame(TermsX2, data = newdata)
  mfZpred2 <- model.frame(TermsZ2, data = newdata)
  
  formYx1 <- reformulate(attr(delete.response(TermsX1), "term.labels"))
  formYz1 <- reformulate(attr(delete.response(TermsZ1), "term.labels"))
  formYx2 <- reformulate(attr(delete.response(TermsX2), "term.labels"))
  formYz2 <- reformulate(attr(delete.response(TermsZ2), "term.labels"))
  X1 <- model.matrix(formYx1, mfXpred1)
  Z1 <- model.matrix(formYz1, mfZpred1)
  X2 <- model.matrix(formYx2, mfXpred2)
  Z2 <- model.matrix(formYz2, mfZpred2)
  ncZ <- ncol(Z1)
  ncZ2 <- ncol(Z2)
  ncX <- ncol(X1)
  
  
  ########## Last time and time to predict
  last.time <- tapply(newdata[[timeVar]], newdata$id, tail, n = 1L)
  n.tp <- sum(!is.na(last.time))
  
  
  
  if (is.null(timesPred)) {
    times.to.pred <- lapply(last.time[!is.na(last.time)], 
                            function (t) seq(t, max(last.time[!is.na(last.time)]) + 1, length = 25L))
  } else {
    times.to.pred <- timesPred
  }
  
  
  ########## Design matrices for prediction of newdata
  if (is.null(timesPred)) {
    data.id2 <- tail(newdata, n = 1L)
    data.id2 <- data.id2[rep(1:nrow(data.id2), 
                             length(times.to.pred)), ]
    data.id2[[timeVar]] <- unlist(times.to.pred)
    ideyes <- as.vector(data.id2$id)
  } else {
    data.id2 <- tail(newdata, n = 1L)
    data.id2 <- data.id2[rep(1:nrow(data.id2), 
                             length(times.to.pred)), ]
    data.id2[[timeVar]] <- unlist(times.to.pred)
    ind2 <- 1:length(times.to.pred)
    data.id2 <- data.id2[ind2,]
    ideyes <- as.vector(data.id2$id)
  }
  
  #data.id2$complication <- getmode(X1[,4])
  
  mfXpred1 <- model.frame(TermsX1, data = data.id2)
  mfZpred1 <- model.frame(TermsZ1, data = data.id2)
  mfXpred2 <- model.frame(TermsX2, data = data.id2)
  mfZpred2 <- model.frame(TermsZ2, data = data.id2)
  
  
  formYx1 <- reformulate(attr(delete.response(TermsX1), "term.labels"))
  formYz1 <- reformulate(attr(delete.response(TermsZ1), "term.labels"))
  formYx2 <- reformulate(attr(delete.response(TermsX2), "term.labels"))
  formYz2 <- reformulate(attr(delete.response(TermsZ2), "term.labels"))
  
  Xpred1 <- model.matrix(formYx1, mfXpred1)
  Zpred1 <- model.matrix(formYz1, mfZpred1)
  Xpred2 <- model.matrix(formYx2, mfXpred2)
  Zpred2 <- model.matrix(formYz2, mfZpred2)
  
  
  ########## Calculate the Empirical Bayes estimates and their (scaled) variance
  betas_y <- apply(sims.list[[1]][,1:ncol(X1)],2,mean) #lmeObject$postMeans$betas2 #fixef(lmeObject)
  betas_x <- apply(sims.list[[2]][,1:ncol(X2)],2,mean) 
  
  sigma_y <- 1/sqrt(mean(sims.list[[3]])) #lmeObject$postMeans$sigma2
  alphas_y <- mean(sims.list[[1]][,(ncol(X1) + 1)])
  D <- solve(matrix(apply(sims.list[[4]],2,mean), ncol(Z1) + ncol(Z2), ncol(Z1) + ncol(Z2)))
  
  
  #id <- match(newdata[[idVar]], unique(newdata[[idVar]]))
  modes <- matrix(0, n.tp, ncol(Z1) + ncol(Z2))
  post_vars <- DZtVinv <- vector("list", n.tp)
  indd <- table(newdata$id)
  
  set.seed(seed)
  
  for (v in seq_len(n.tp)) {
    
    scale = 1.6
    
    y_y <- y1 #model.response(mfX_y)
    
    y_x <- y2 #model.response(mfX_x)
    
    betas.new_y <- betas_y
    betas.new_x <- betas_x
    sigma.new <- sigma_y
    D.new <- D
    alpha.new <- alphas_y
    
    
    ff <- function(b, y_y, y_x, X_new_y, Z_new_y, X_new_x, Z_new_x, betas.new_y, betas.new_x, sigma.new, 
                   D.new, 
                   alpha.new) {
      - ( -(1/(2*sigma.new^2))*(t(y_y-X_new_y%*%betas.new_y-Z_new_y%*%b[(1:ncZ)] - alpha.new * (X_new_x%*%betas.new_x + Z_new_x%*% b[(ncZ+1):(ncZ+ncZ2)])   )%*%(y_y-X_new_y%*%betas.new_y-Z_new_y%*%b[(1:ncZ)] -alpha.new * (X_new_x%*%betas.new_x + Z_new_x%*% b[(ncZ+1):(ncZ+ncZ2)])  )  ) +
            sum(like(pii(X_new_x, betas.new_x, Z_new_x, b[(ncZ+1):(ncZ+ncZ2)]), y_x)) -
            (1/2*(t(b)%*%solve(D.new)%*%(b)))  )  }
    opt <- try(optim(rep(0, (ncol(Z1) + ncol(Z2))), ff, y_y = y_y[], y_x = y_x[], 
                     X_new_y = X1[,], Z_new_y = Z1[,], X_new_x = X2[,], Z_new_x = Z2[,],
                     betas.new_y = betas.new_y, betas.new_x = betas.new_x, sigma.new = sigma.new, 
                     D.new = D.new,
                     alpha.new = alpha.new, 
                     method = "BFGS", hessian = TRUE), TRUE)
    
    if (inherits(opt, "try-error")) {
      gg <- function(b, y_y, y_x, X_new_y, Z_new_y, X_new_x, Z_new_x, betas.new_y, betas.new_x, sigma.new, 
                     D.new, 
                     alpha.new) {cd(b, ff, y_y = y_y, y_x = y_x,  
                                    X_new_y = X1[,], Z_new_y = Z1[,], X_new_x = X2[,], Z_new_x = Z2[,],
                                    betas.new_y = betas.new_y, betas.new_x = betas.new_x, sigma.new = sigma.new, 
                                    D.new = D.new,
                                    alpha.new = alpha.new)}
      opt <- optim(rep(0, (ncol(Z1) + ncol(Z2))), ff, gg, y_y = y_y, y_x = y_x, 
                   X_new_y = X1[,], Z_new_y = Z1[,], X_new_x = X2[,], Z_new_x = Z2[,],
                   betas.new_y = betas.new_y, betas.new_x = betas.new_x, sigma.new = sigma.new, 
                   D.new = D.new,
                   alpha.new = alpha.new, 
                   method = "BFGS", hessian = TRUE)
    }
    modes[v, ] <- opt$par
    DZtVinv[[v]] <- opt$hessian/scale
    post_vars[[v]] <- scale * solve(opt$hessian)        
  }
  
  
  

  
  ########## Predictions
  ypred2 <- as.vector(c(Xpred2 %*% betas_x) + rowSums(Zpred2 * matrix(rep(modes[(ncZ + 1):(ncZ + ncZ2)], each = nrow(Zpred2)), ,ncol(Zpred2))))
  
  y22 <- rbinom(length(ypred2), 1, exp(ypred2)/(1+exp(ypred2)))
  ypred1 <- as.vector(c(Xpred1 %*% betas_y) + rowSums(Zpred1 * matrix(rep(modes[(1):(ncZ)], each = nrow(Zpred1)), ,ncol(Zpred1))) +
                        alphas_y * y22)
  
  ########## Output
  out <- data.frame(ypred1 = unlist(ypred1), 
                    ypredIN = unlist(ypred2), 
                    times.to.pred = unlist(times.to.pred), ideye = unlist(ideyes))
  out
}



