
#### functions 
blasso <- function(b, lam) {
  max(b - lam, 0)
}
bqHALF <- function(b, lam, upper = T, eps = 0) {
  bOLSmin <- (lam/4)^(2/3) + lam/2 * ((lam/4)^(2/3))^(-1/2)
  bsign <- sign(b)
  lam <- complex(real = lam)
  I <- complex(real = 0, imaginary = 1)
  
  #if ((b < bOLSmin) & (eps != 0)) {
  #  b <- complex(real = abs(b), imaginary = eps)
  #} else {
  #  if (b < bOLSmin) {
  #    return (0)
  #  }
  #  b <- complex(real = abs(b), imaginary = eps)
  #}
  if (b < 0.998*bOLSmin) {
    return (NaN)
  }
  
  term1 <-  2 * b / 3
  term2 <- -(1 - I * sqrt(3)) * b^2 / (3 * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3))
  term3 <- -1/12 * (1 + I * sqrt(3)) * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3)
  
  if (!upper) {
    term2 <- -(1 + I * sqrt(3)) * b^2 / (3 * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3))
    term3 <- -1/12 * (1 - I * sqrt(3)) * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3)
  }
  
  x <- term1 + term2 + term3
  
  if (Re(x) < 0) {
    return (0)
  } else {
    return (bsign * x)
  }
  #if (abs(Im(x)) > 1e-3) {
  #  return (0)
  #} else {
  #  return (Re(bsign * x))
  #}
}
bqHALF_shift <- function(b, lam, upper = T, eps = 0.01) {
  bsign <- sign(b)
  lam <- complex(real = lam)
  I <- complex(real = 0, imaginary = 1)
  
  b <- complex(real = abs(b), imaginary = eps)
  term1 <-  2 * b / 3
  term2 <- -(1 - I * sqrt(3)) * b^2 / (3 * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3))
  term3 <- -1/12 * (1 + I * sqrt(3)) * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3)

  if (upper) { # we actually do the opposite to get the upper branch
    term2 <- -(1 + I * sqrt(3)) * b^2 / (3 * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3))
    term3 <- -1/12 * (1 - I * sqrt(3)) * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3)
  }
  
  
  x <- term1 + term2 + term3
  
  if (Re(x) < 0) {
    return (0)
  } else {
    return (bsign * x)
  }
}

eps <- 0.1
bmax <- 2
lam <- 1

beta_ols <- seq(0, bmax, length.out = 1e3)
beta_ols2 <- beta_ols # hacky solution to plot x ~ x
beta_q_up <- sapply(beta_ols, function(b) bqHALF(b = b, lam = lam, upper = T, eps = 0))
beta_q_lo <- sapply(beta_ols, function(b) bqHALF(b = b, lam = lam, upper = F, eps = 0))
beta_q_s1 <- sapply(beta_ols, function(b) bqHALF_shift(b = b, lam = lam, upper = T, eps = eps))
beta_q_s2 <- sapply(beta_ols, function(b) bqHALF_shift(b = b, lam = lam, upper = F, eps = eps))
beta_lasso <- sapply(beta_ols, function(b) blasso(b, lam))

plot(NA, 
     xlab = expression(beta),
     ylab = expression(beta^OLS),
     xlim = c(min(beta_ols), max(beta_ols)), 
     ylim = c(min(Re(beta_q_up), beta_lasso), max(Re(beta_q_up), beta_lasso)),
     main = substitute(
       paste("q = 0.5, ", lambda, " = ", lam, ", ", epsilon, " = ", eps), 
       list(lam = lam, eps = eps)))
lines(beta_ols ~ beta_ols2, lwd = 1, col = 'gray70')
lines(beta_lasso ~ beta_ols, lwd = 1, col = 'red', lty = 'dashed')
lines(Re(beta_q_up) ~ beta_ols, lwd = 2)
lines(Re(beta_q_lo) ~ beta_ols, lwd = 2)
lines(Re(beta_q_s1) ~ beta_ols, lwd = 1.5, lty = 'dotted', col = 'green3')
lines(Re(beta_q_s2) ~ beta_ols, lwd = 1.5, lty = 'dotdash', col = 'green3')
legend("topleft", col = c("black", "green3", "green3", "red", "gray70"),
       lwd = c(2, 1.5, 1.5, 1, 1), bty = 'n', seg.len = 2,
       lty = c("solid", "dotted", "dotdash", "dashed", "solid"),
       legend = c("1/2-norm no shift", "1/2-norm pos shift", "1/2-norm neg shift", "Lasso", "OLS"))


plot(NA, 
     xlab = expression(beta),
     ylab = expression(beta^OLS),
     xlim = c(min(beta_ols), max(beta_ols)), 
     ylim = c(min(Re(beta_q), beta_lasso), max(Re(beta_q), beta_lasso)),
     main = substitute(
       paste("q = 0.5, ", lambda, " = ", lam, ", ", epsilon, " = ", eps), 
       list(lam = lam, eps = eps)))
lines(beta_ols ~ beta_ols2, lwd = 1, col = 'gray70')
lines(beta_lasso ~ beta_ols, lwd = 1, col = 'red', lty = 'dashed')
lines(Re(beta_q_s1) ~ beta_ols, lwd = 2, lty = 'dotted', col = 'green3')
lines(Re(beta_q_s2) ~ beta_ols, lwd = 2, lty = 'dotdash', col = 'green3')
legend("topleft", col = c("green3", "green3", "red", "gray70"),
       lwd = c(2, 2, 1, 1), bty = 'n', seg.len = 2,
       lty = c("dotted", "dotdash", "dashed", "solid"),
       legend = c("1/2-norm pos shift", "1/2-norm neg shift", "Lasso", "OLS"))



