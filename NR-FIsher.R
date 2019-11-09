Fisher.it <-  function(Y, X, pi0, niter=1, print=F) {
  pi <-  pi0
  for (i in 1:niter) {
    W <-  pi*(1-pi)
    Z <- log(pi/(1-pi)) + (Y - pi)/(pi*(1-pi))
    lmobj <-  lm(Z ~ X - 1, weights=W)
    beta <-  lmobj$coef
    eta <-  X %*% beta
    pi <-  exp(eta)/(1 + exp(eta))
    if (print) {
      print(paste("Iteration ", as.character(i), ": Betahat"))
      print(beta)}
  }
  XWX <- t(lmobj$R) %*% lmobj$R
  return(beta, XWX, pi, W)
}

