X<- matrix(rnorm(20), ncol = 2)
y <- sample(c(0,1), replace=TRUE, size=10)

m <- rep(1, length(y))
X <- cbind(X, c(1,1,1,1,1,1,1,1,1,1))
r <- ncol(X) - 1 # number of regression coefficients - 1
beta_new <- c(log(sum(y) / sum(m - y)), rep(0, r))

# Define logistic regression
logistic <- function(X, beta) {
  p <- exp(X %*% beta)/(1+exp(X %*% beta))
}

Fisher_scoring <- function(X, y, m, beta_new, maxiter, delta.betas) {
  # X - a matrix of covariates, must contain a column of 1s as an intercept
  # y - a vector of "successes" 
  # m - a vector denoting the total number of trials - equals to 1 for Bernoulli.
  #  m is equal to 1 as in 2a). 
  # beta_new - a vector of initial guess for betas
  # maxiter - a bound for a maximum number of iterations
  # delta.beta - the bound for difference between the consecutive 
  approximated betas that indicates convergence
 
  beta_old <- rep(Inf, length(beta_new))
  # Euclidian distance 
  diff.betas <- sqrt(sum((beta_new - beta_old)^2))
  iter <- 1 
  
  while((iter<=maxiter) & (diff.betas>delta.betas)) {
    iter <- iter + 1
    beta_old <- beta_new
    # estimate the improvement part of the alogrithm: (t(X)WX)^-1*t(X)(y-mp)
    W <- diag(as.vector(m*logistic(X, beta_old)*(1-logistic(X, beta_new))))
    mp <- m * logistic(X, beta_old)
    # XWX is calculated with Tikhonov reguliarisation -- multicolinearity problem
    XWX <- t(X) %*% W %*% X + diag(1, 58, 58) * 0.0001
    Fscore <- t(X) %*% (y-mp)
    # (XWX)^-1 * Fscore = Z, in R - solve XWX * Z = Fscore
    improvement <- solve(XWX, Fscore)
    # the final algorithm of iterations
    beta_new <- beta_old + improvement
    
    diff.betas <- sqrt(sum((beta_new - beta_old)^2))
  }
  
  result <- list()
  result$beta.MLE <- beta_new
  result$iterations <- iter - 1
  return(result)
}

set.seed(1620789)
Fisher_scoring(X, y, m, beta_new, maxiter=50, delta.beta = 0.0001)

# Comparison with the inbuilt model
glm.fit <- glm( y~X, family = binomial)
summary(glm.fit)
# same result, slower convergence.

