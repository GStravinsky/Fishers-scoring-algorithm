######################################
# Newton Raphson and Fisher scoring #
#####################################

# Define logistic regression
logistic <- function(X, beta) {
  p <- exp(X %*% beta)/(1+exp(X%*%beta))
}

Fisher_scoring <- function(X, y, m, beta1, iter, delta.beta) {
  
  # setting ground for the algirithm
  # beta2 is used instrumentally to generate comparisons and calculate differences
  beta2 <- rep(Inf, length(beta1))
  # Euclidian distance to calulate the improvements
  diff.betas <- sqrt(sum(beta1-beta2)^2)
  i <- 1 # first interation
  # store there the approximations
  beta.hist <- matrix(beta1, nrow=1)
  
  # Actual process begins here
  while((i<=iter) & (diff.betas>delta.beta)) {
    # keeping count of iterations
    i <- i + 1
    # because the original input beta1 is the initial guess, 
    # and beta1 in the code bellow represents improved the approximation of beta 
    beta2 <- beta1
    # estimate the improvement part of the alogrithm: (t(X)WX)^-1*t(X)(y-mp)
    W <- diag(as.vector(m*logistic(X, beta2)*(1-logistic(X, beta2))))
    mp <- m * logistic(X, beta2)
    # (t(X)WX)^-1*t(X)(y-mp) = Z, in R - solve t(X)WX * Z = t(X)(y-mp)
    improvement <- solve(t(X) %*% W %*% X, t(X) %*% (y-mp))
    # the final algorithm of iterations
    beta1 <- beta2 + improvement
    
    diff.betas <- sqrt(sum(beta1-beta2)^2)
    beta.hist <- rbind(beta.hist, (matrix(beta1, nrow=1)))
  }
  
  result <- list()
  result$beta.hat <- beta1
  result$beta.hist <- beta.hist
  result$iter <- i - 1
  return(result)
}

##########################
# DUMMY DATA FOR A TRIAL #
##########################


## Beetles data set
# conc = CS2 concentration
# y = number of beetles killed
# n = number of beetles exposed
# rep = Replicate number (1 or 2)
beet <- read.table("http://statacumen.com/teach/SC1/SC1_11_beetles.dat", header = TRUE)
beet$rep <- factor(beet$rep)
# create data variables: m, y, X
n <- nrow(beet)
m <- beet$n
y <- beet$y
X.temp <- beet$conc
# quadratic model
X <- matrix(c(rep(1,n), X.temp, X.temp^2), nrow = n)
colnames(X) <- c("Int", "conc", "conc2")
r <- ncol(X) - 1 # number of regression coefficients - 1

# initial beta vector
beta1 <- c(log(sum(y) / sum(m - y)), rep(0, r))

Fisher_scoring(X, y, m, beta1, iter=10, delta.beta=0.001 )
