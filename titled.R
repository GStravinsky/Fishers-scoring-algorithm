#Titled.R By GS
#TODO: Define logistic model
#TODO: Create x and b vectors
#TODO: define log likelihood
#TODO: FIRST derivative
#TODO: Second derivative - Fishers information
#TODO: Plug it all into the algorithm

logistic <- function(X, beta) {
  # compute vector p of probabilities for logistic regression with logit link
  # X - covaroates, beta - coefficients, p - probability of being in category 1
  X <- as.matrix(X)
  beta <- as.vector(beta)
  p <- exp(X %*% beta) / (1 + exp(X %*% beta))
  return(p)
}

loglike_logistic <-  function(y, p) {
  # binomial log likelihood function
  # input: vectors: y = groups; p = probabilities
  # output: log-likelihood l, a scalar
  l <- t(y) %*% log(p) + t(1-y) %*% log(1 - p)
  return(l)
}
# why do you need to transpose y???


Newton.Raphson <- function(X, y, beta1, eps1=0.01, eps2=0.01, maxit=50) {
  #X := nX(r+1) explanatory variable matrix
  # y := column vector of "sucess" counts
  #  beta.1 = (r+1)-by-1 vector of starting values 
  # r - number of covariates
  # eps1 = absolute convergence criterion for beta
  # eps2 = absolute convergence criterion for log-likelihood
  # maxit = maximum allowable number of iterations
  
  beta.2 <- rep(-Inf, length(beta.1)) # init beta.2
  diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
  
  #Construct likelihood depending on you data X, and chosen initial beta1
  llike.1 <- loglike_logistic(y, logistic(X, beta.1))
  llike.2 <- loglike_logistic(y, logistic(X, beta.2))
  diff.like <- abs(llike.1 - llike.2) # difference between likelihoods
  if (is.nan(diff.like)) { diff.like <- 1e9 }
  
  i <- 1 # initial iteration index
  
  
  NR.hist <- data.frame(i, diff.beta, diff.like, llike.1, step.size = 1) # iteration history
  beta.hist <- matrix(beta.1, nrow = 1)
  # they measure distance in euclidian space. could make sense for a multivariate case. 
  
  # whilst there are fewer than max iterations and there is some improvement in beta (getting closer to zero) as well as likelihood estimation (getting bigger)
  while ((i <= maxit) & (diff.beta > eps1) & (diff.like > eps2)) {
    i <- i + 1 # increment iteration
    
    # update beta
    beta.2 <- beta.1 # old guess is current guess
    
    # The mor complicated components of the interation  
    W <- diag(as.vector(logistic(X, beta.2)*(1-logistic(X, beta.2))))
    # Score - first derivative.
    score <- t(X) %*% (y-logistic(X, beta.2))
    # Fisher information - inverse
    FI <- solve(t(X) %*% v.2 %*% X)
    improvement <- FI %*% score
    
    beta.1 <- beta.2 - improvement
    
  }
  
}
