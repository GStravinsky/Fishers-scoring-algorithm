
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

beta.1 <- c(log(sum(y) / sum(m - y)), rep(0, r))

logistic <- function(X, beta) {
  # compute vector p of probabilities for logistic regression with logit link
  # X - covaroates, beta - coefficients, p - probability of being in category 1
  X <- as.matrix(X)
  beta <- as.vector(beta)
  p <- exp(X %*% beta) / (1 + exp(X %*% beta))
  return(p)
}

loglike_logistic <-  function(y, m, p) {
  # binomial log likelihood function
  # input: vectors: y = groups; p = probabilities
  # output: log-likelihood l, a scalar
  l <- t(y) %*% log(p) + t(m-y) %*% log(1 - p)
  return(l)
}
# why do you need to transpose y???


one.iter <- function(X, y, m, beta.1){
# update beta
beta.2 <- beta.1 # old guess is current guess
mu.2 <- m * logistic(X, beta.2)
# The mor complicated components of the interation  
W <- diag(as.vector(m*logistic(X, beta.2)*(1-logistic(X, beta.2))))
# Score - first derivative.
score <- t(X) %*% (y-mu.2)
# Fisher information - inverse
FI <- solve(t(X) %*% W %*% X)
improvement <- FI %*% score

beta.1 <- beta.2 - improvement
return(beta.1)
}

out <- one.iter(X, y, m, beta.1)

