
#### set up
library(genlasso) # genlass()
library(MASS) # lm.ridge()

#### functions 
bqHALF <- function(b, lam, upper = T) {
  bsign <- sign(b)
  b <- complex(real = abs(b))
  lam <- complex(real = lam)
  I <- complex(real = 0, imaginary = 1)
  
  term1 <-  2 * b / 3
  term2 <- -(1 - I * sqrt(3)) * b^2 / (3 * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3))
  term3 <- -1/12 * (1 + I * sqrt(3)) * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3)
  
  if (!upper) {
    term2 <- -(1 + I * sqrt(3)) * b^2 / (3 * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3))
    term3 <- -1/12 * (1 - I * sqrt(3)) * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3)
  }
  
  x <- term1 + term2 + term3
  
  if (abs(Im(x)) > 1e-3) {
    return (0)
  } else {
    return (Re(bsign * x))
  }
}

#### generate data
set.seed(125)
correlated <- F; sig <- 0
nsims <- 2e3
err <- 1e-3 # if coef < err map it to zero
p <- 5
n <- 1e2
lams <- 10^seq(-4, 1, length.out = 5e3)

betas <- rnorm(p)
D <- diag(p)

pt <- proc.time()
sim <- replicate(nsims, {
  Z <- matrix(rnorm(n * p), ncol = p)
  eps <- rnorm(n)

  #X <- qr.Q(qr(Z)) # generate orthgonal data
  # turns out this is 'close enough' to orthogonal (when observations are iid)
  Zstd <- apply(Z, 2, function(z) 1/sqrt(n - 1) * (z - mean(z))/sd(z)) 
  Xstd <- Zstd

  ## generate correlated normals
  if (correlated) {
    C <- matrix(0, ncol = p, nrow = p)
    diag(C) <- rnorm(p, mean = 1, sd = sig); C[lower.tri(C)] <- rnorm(p * (p - 1) / 2, sd = sig)
    SIG <- C %*% t(C) # covariance matrix
    RHO <- (1/sqrt(diag(SIG)) * diag(p)) %*% SIG %*% (1/sqrt(diag(SIG)) * diag(p)) # correl mat
    
    Xtmp <- t(C %*% t(Z))
    Xstd <- apply(Xtmp, 2, function(x) 1/sqrt(n - 1) * (x - mean(x))/sd(x))
  }

  y <- Xstd %*% betas + eps
  
  #### models
  # ols estimates
  bols <- coef(lm(y ~ Xstd))[-1] # exclude intercept
  bols <- bols[!is.na(bols)]
  
  # lasso estimates
  cf_lasso <- coef(genlasso(y, Xstd, D), lambda = lams)
  beta_lasso <- apply(cf_lasso$beta, 2, rev)
  
  # ridge estimates
  #beta_ridge <- coef(lm.ridge(y ~ Xstd, lambda = lams))[,-1] # rid of intercept
  
  # get 1/2-norm estimates
  bq_ests_list <- lapply(lam_lass1, FUN = function(l) sapply(bols, FUN = function(b) bqHALF(b, l, upper = T)))
  beta_q <- matrix(unlist(bq_ests_list), ncol = length(bols), byrow = T)
  
  lam_lasso_zero <- lams[apply(beta_lasso, 1, function(b) which(abs(b) < err)[1])]
  lam_q_zero <- lams[apply(beta_q, 2, function(b) which(abs(b) < err)[1])]
  
  rbind(bols, lam_lasso_zero, lam_q_zero)
})
proc.time() - pt

bols_sim <- sim[1,,] # bols
lam_lasso_sim <- sim[2,,] # lam_lasso_zero
lam_q_sim <- sim[3,,] # lam_q_zero

k <- 2
hist(bols_sim[k,], freq = F)
abline(v = betas[k], lwd = 2, col = 'red')
lines(density(bols_sim[k,]), lwd = 2)

hist(lam_lasso_sim[k,], freq = F)
lines(density(lam_lasso_sim[k,]), lwd = 2)

hist(lam_q_sim[k,], freq = F)
lines(density(lam_q_sim[k,]), lwd = 2)

plot(lam_q_sim[k,] ~ lam_lasso_sim[k,])










