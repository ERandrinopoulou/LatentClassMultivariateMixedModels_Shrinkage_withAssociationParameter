model
{
  for (i in 1:N) {
    for (j in offset[i]:(offset[i + 1] - 1)) {
      muy[j] <- inprod(betas[1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncZ], Z[j, 1:ncZ]) + betas[ncX + 1] * y2[j] + betas[ncX + 2] * y3[j]
      y[j] ~ dnorm(muy[j], tau)
    }
    for (j in offset[i]:(offset[i + 1] - 1)) {
      muy2[j] <- inprod(betas2[1:ncX2], X2[j, 1:ncX2]) + inprod(b[i, (ncZ + 1):(ncZ + ncZ2)], Z2[j, 1:ncZ2])
      Pr2[j] <- max(1.00000E-05, min(0.99999, (exp(muy2[j])/(1 + exp(muy2[j])))))
      y2[j] ~ dbin(Pr2[j], 1)
    }
    for (j in offset[i]:(offset[i + 1] - 1)) {
      muy3[j] <- inprod(betas3[1:ncX3], X3[j, 1:ncX2]) + inprod(b[i, (ncZ + ncZ2 + 1):(ncZ + ncZ2 + ncZ3)], Z3[j, 1:ncZ3])
      Pr3[j] <- max(1.00000E-05, min(0.99999, (exp(muy3[j])/(1 + exp(muy3[j])))))
      y3[j] ~ dbin(Pr3[j], 1)
    }
    b[i, 1:nb] ~ dmnorm(mu0[], inv.D[, ])
  }

    # Prior for betas
    for(j in 1:(ncX + 2)){
        betas[j] ~ dnorm(0,inv.var.b1)
    }
    for(j in 1:(ncX2)){
        betas2[j] ~ dnorm(0,inv.var.b2)
    }
    for(j in 1:(ncX3)){
        betas3[j] ~ dnorm(0,inv.var.b3)
    }

  
    inv.var.b1 ~ dgamma(ni/2, 1/((ni/2) * sig^2))
    inv.var.b2 ~ dgamma(ni/2, 1/((ni/2) * sig^2))
    inv.var.b3 ~ dgamma(ni/2, 1/((ni/2) * sig^2))
    niInv ~ dunif(0.00000E+00, 1)
    ni <- 1/niInv
    sig ~ dunif(0.00000E+00, 100)

  tau ~ dgamma(priorA.tau, priorB.tau)
  inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)
}
