logistic <- function(X, beta) {
  # compute vector p of probabilities for logistic regression with logit link
  # X - covaroates, beta - coefficients, p - probability of being in category 1
  X <- as.matrix(X)
  beta <- as.vector(beta)
  p <- exp(X %*% beta) / (1 + exp(X %*% beta))
  return(p)
}

loglike_logistic <-  function(y, p, m) {
  # binomial log likelihood function
  # input: vectors: y = groups; p = probabilities
  # output: log-likelihood l, a scalar
  l <- t(y) %*% log(p) + t(1-y) %*% log(1 - p)
  return(l)
}

Fisher <- function(Y, X, beta.1, m, niter){
  for(i in 1:niter){
    beta.2 <- beta.1
    W <- diag(as.vector(m*logistic(X, beta.2)*(1-logistic(X, beta.2))))# Score - first derivative.
    mu.2 <- m * logistic(X, beta.2) # m * p is mean
    score <- t(X) %*% (y-mu.2)
    # Fisher information - inverse
   # XWX <- solve(t(X) %*% W %*% X)
    improvement <- solve(t(X) %*% W %*% X, score) # solve for increment
    beta.1 <- beta.2 - improvement
  }
  return(beta.1)
}

Fisher(y, X, beta.1, m, 10)
