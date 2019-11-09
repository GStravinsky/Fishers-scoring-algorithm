X<- matrix(rnorm(20), ncol = 2)
y <- sample(c(0,1), replace=TRUE, size=10)

m <- rep(1, length(y))
X <- cbind(X, c(1,1,1,1,1,1,1,1,1,1))
r <- ncol(X) - 1 # number of regression coefficients - 1
beta_new <- c(log(sum(y) / sum(m - y)), rep(0, r))


# Define logistic regression
logistic <- function(X, beta) {
  p <- exp(X %*% beta)/(1+exp(X%*%beta))
}

Fisher_scoring <- function(X, y, m, beta_new, maxiter, delta.beta) {
  
  # beta2 is used instrumentally to generate comparisons and calculate differences
  beta_old <- rep(Inf, length(beta_new))
  # Euclidian distance to calulate the improvements
  diff.betas <- sqrt(sum(beta_new-beta_old)^2)
  iter <- 1 
  beta.hist <- matrix(beta_new, nrow=1)
  
  while((iter<=iter) & (diff.betas>delta.beta)) {
    iter <- iter + 1
    # because the original input beta1 is the initial guess, 
    # and beta1 in the code bellow represents improved the approximation of beta 
    beta_old <- beta_new
    # estimate the improvement part of the alogrithm: (t(X)WX)^-1*t(X)(y-mp)
    W <- diag(as.vector(m*logistic(X, beta_old)*(1-logistic(X, beta_new))))
    mp <- m * logistic(X, beta_old)
    # (t(X)WX)^-1*t(X)(y-mp) = Z, in R - solve t(X)WX * Z = t(X)(y-mp)
    improvement <- solve(t(X) %*% W %*% X, t(X) %*% (y-mp))
    # the final algorithm of iterations
    beta_new <- beta_old + improvement
    
    diff.betas <- sqrt(sum(beta_new-beta_old)^2)
    beta.hist <- rbind(beta.hist, (matrix(beta_new, nrow=1)))
  }
  
  result <- list()
  result$beta.MLE <- beta_new
  result$beta.history <- beta.hist
  result$iterations <- iter - 1
  return(result)
}

Fisher_scoring(X, y, m, beta_new, maxiter=50, delta.beta = 0.0001)

glm.fit <- glm( y~X, family = binomial)
summary(glm.fit)

