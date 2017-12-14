
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
exp_xaxis <- T
p <- 10
n <- 1e3
Z <- matrix(rnorm(n * p), ncol = p)
betas <- rnorm(p)
eps <- rnorm(n)

#X <- qr.Q(qr(Z)) # generate orthgonal data
# turns out this is 'close enough' to orthogonal (when observations are iid)
Zstd <- apply(Z, 2, function(z) 1/sqrt(n - 1) * (z - mean(z))/sd(z)) 
Xstd <- Zstd

## generate correlated normals
C <- matrix(0, ncol = p, nrow = p)
diag(C) <- rnorm(p, mean = 1, sd = sig); C[lower.tri(C)] <- rnorm(p * (p - 1) / 2, sd = sig)
SIG <- C %*% t(C) # covariance matrix
RHO <- (1/sqrt(diag(SIG)) * diag(p)) %*% SIG %*% (1/sqrt(diag(SIG)) * diag(p)) # correl mat
if (correlated) {
  Xtmp <- t(C %*% t(Z))
  Xstd <- apply(Xtmp, 2, function(x) 1/sqrt(n - 1) * (x - mean(x))/sd(x))
}

y <- Xstd %*% betas + eps
D <- diag(p)

round(t(Xstd) %*% Xstd, 3)[1:min(p, 10), 1:min(p, 10)]
round(RHO, 3)
round(cor(Xstd), 3)

# ols estimates
bols <- coef(lm(y ~ Xstd))[-1] # exclude intercept
bols <- bols[!is.na(bols)]

betas
bols

#### models
lams <- 10^seq(-3, 4, length.out = 1e3)
lams <- 10^log10(c(seq(0.001, 0.999, 0.01), seq(1, 10, 1), seq(15, 50, 5), seq(60, 1e4, 10)))

# lasso estimates
mod_lass1 <- genlasso(y, Xstd, D)
cf_lass1 <- coef(mod_lass1, lambda = lams)

beta_lass1 <- apply(cf_lass1$beta, 2, rev)
lam_lass1 <- as.numeric(colnames(cf_lass1$beta))

# ridge estimates
beta_ridge <- coef(lm.ridge(y ~ Xstd, lambda = lam_lass1))[,-1] # rid of intercept

# get 1/2-norm estimates
bq_ests_list <- lapply(lam_lass1, FUN = function(l) sapply(bols, FUN = function(b) bqHALF(b, l, upper = T)))
beta_q <- matrix(unlist(bq_ests_list), ncol = length(bols), byrow = T)

#### plots
lam_lass1_map <- -1 * (2/pi * atan(lam_lass1) - 1) # map lambda to [0, 1] for visualization purposes
#lam_lass1_map <- 1 - tanh(lam_lass1)
#lam_lass1_map <- 2 - 2 * pnorm(lam_lass1) 
#lam_lass1_map <- 2 - 2 * pcauchy(lam_lass1)
#lam_lass1_map <- 2 - 2/(1 + exp(-lam_lass1))

### plot all together
ylm <- c(min(beta_lass1, beta_q), max(beta_lass1, beta_q))
#ylm <- c(-2, 10)
xlm <- c(0, 1.02)
ols_pos <- 1
if (exp_xaxis) {
  lam_lass1_map <- log10(lam_lass1)
  xlm <- rev(c(min(lam_lass1_map[is.finite(lam_lass1_map)]), 
               max(lam_lass1_map[is.finite(lam_lass1_map)]) * 1.02))
  ols_pos <- min(xlm)
}
plot(beta_lass1[1,] ~ lam_lass1_map, type = 'l', lwd = 2, col = 'darkblue',
     ylim = ylm, las = 1, lty = 'longdash',
     xlim = xlm,
     xlab = expression(-(2/pi * atan(lambda) - 1)), 
     ylab = expression(beta),
     main = "Lasso vs. 1/2-norm\nCoefficient Estimates")
for (i in 2:nrow(beta_lass1)) {
  lines(beta_lass1[i,] ~ lam_lass1_map, lwd = 2, col = 'darkblue', lty = 'longdash')
}
for (i in 1:ncol(beta_ridge)) {
  lines(beta_ridge[,i] ~ lam_lass1_map, lwd = 2, col = "green3", lty = 'dotted')
}
for (i in 1:ncol(beta_q)) {
  lines(beta_q[,i] ~ lam_lass1_map, lwd = 2, col = 'red')
  text(x = 1.03 * ols_pos, y = beta_q[1,i], labels = as.character(i))
}
abline(v = ols_pos, col = 'gray70')
points(x = rep(ols_pos, length(bols)), y = bols, pch = 19, col = "gray70")
legend("topleft", legend = c("Lasso", "Ridge", "1/2-norm", "OLS"), 
       col = c("darkblue", "green3", "red", "gray70"), lty = c("longdash", "dotdash", "solid", "solid"), 
       lwd = c(2, 2, 2, 1), pch = c(NA, NA, NA, 19), bty = 'n', seg.len = 2.5)

############## TESTS ###############
### just lasso plots
# plot(beta_lass1[1,] ~ lam_lass1_map, type = 'l', lwd = 2, col = 'darkblue',
#      ylim = c(min(beta_lass1), max(beta_lass1)), las = 1,
#      xlim = c(0, 1.04),
#      xlab = expression(-(2/pi * atan(lambda) - 1)), ylab = "", main = "Lasso Coefficient Estimates")
# text(x = min(lam_lass1_map), y = (max(beta_lass1) + min(beta_lass1))/2, adj = 8, 
#      labels = expression(beta), xpd = TRUE)
# text(x = 1.04, y = beta_q[1,1], labels = as.character(i))
# for (i in 2:nrow(beta_lass1)) {
#   lines(beta_lass1[i,] ~ lam_lass1_map, lwd = 2, col = 'darkblue')
#   text(x = 1.04, y = beta_q[1,i], labels = as.character(i))
# }
# points(x = rep(1, length(bols)), y = bols, pch = 19, col = "gray50")
# legend("topleft", legend = c("Lasso", "1/2-norm", "OLS"), 
#        col = c("darkblue", "red", "gray50"), lty = c("longdash", "solid", NA), 
#        lwd = c(2, 2, NA), pch = c(NA, NA, 19), bty = 'n')

