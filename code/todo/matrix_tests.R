set.seed(124)
n <- 1e2
beta_true <- c(1, 5, -3, 0, 0, 1)
sigma <- 1
p <- length(beta_true)

X <- matrix(rnorm(n * p), ncol = p)
Y <- X %*% beta_true + rnorm(n, sd = sigma)

beta_q_est <- function(beta, x, y, q, lambda) {
  svdX <- svd(X)
  U <- svdX$u
  D <- svdX$d
  V <- svdX$v
  
  Dstar <- diag((D/(D^2 + lambda * q / 2 * abs(beta)^(q - 2))))
  
  sum((beta - V %*% Dstar %*% t(U) %*% Y)^2)
}

q <- 1
lambda_v <- seq(20, 150, length.out = 50)
beta_est <- matrix(nrow = length(lambda_v), ncol = p)

for (i in 1:length(lambda_v)) {
  lambda <- lambda_v[i]
  beta_est[i,] <- optim(rep(1, p), function(b) beta_q_est(b, x = X, y = Y, q = q, lambda = lambda), 
                        method = "CG")$par
}

lambda_v_map <- -2/pi * atan(lambda_v) + 1

plot(beta_est[,1] ~ lambda_v_map, type = 'l', ylim = c(min(beta_est), max(beta_est)),
     lwd = 2, col = 'darkblue', log = 'x')
for (j in 2:ncol(beta_est))
  lines(beta_est[,j] ~ lambda_v_map, col = 'darkblue', lwd = 2)


#round(beta_test, 4) - beta_true
#round(beta_test, 4)










