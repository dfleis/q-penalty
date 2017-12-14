#### functions 
blasso <- function(b, lam) {
  max(b - lam, 0)
}
bqHALF <- function(b_in, lam_in, upper = T, eps = 0) {
  bOLSmin <- (lam/4)^(2/3) + lam/2 * ((lam/4)^(2/3))^(-1/2)
  if (b_in < bOLSmin) 
    return (NaN)
    
  bsign <- sign(b_in)
  lam <- complex(real = lam_in)
  I <- complex(real = 0, imaginary = 1)
  b <- complex(real = b_in)
  
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

eps <- 0.1
bmax <- 3
lam <- 1

beta_ols <- seq(0, bmax, length.out = 1e4)
beta_ols2 <- beta_ols # hacky solution to plot x ~ x
beta_q_up <- sapply(beta_ols, function(b) bqHALF(b_in = b, lam_in = lam, upper = T, eps = 0))
beta_q_lo <- sapply(beta_ols, function(b) bqHALF(b_in = b, lam_in = lam, upper = F, eps = 0))
beta_lasso <- sapply(beta_ols, function(b) blasso(b, lam))

plot(NA, 
     xlab = expression(beta),
     ylab = expression(beta^OLS),
     xlim = c(min(beta_ols), max(beta_ols)), 
     ylim = c(0, max(Re(beta_q_up[!is.nan(beta_q_up)]))),
     main = substitute(
       paste("q = 0.5, ", lambda, " = ", lam, ", ", epsilon, " = ", eps), 
       list(lam = lam, eps = eps)))
lines(beta_ols ~ beta_ols2, lwd = 0.75, col = 'gray70')
lines(Re(beta_q_up) ~ beta_ols, lwd = 2, lty = 'solid', col = 'blue4')
lines(Re(beta_q_lo) ~ beta_ols, lwd = 2, lty = 'solid', col = 'deepskyblue')
legend("topleft", col = c("blue4", "deepskyblue", "gray70"),
       lwd = c(2, 2, 0.75), bty = 'n', seg.len = 2,
       lty = c("solid", "solid", "solid"),
       legend = c("1/2-norm top branch", "1/2-norm top branch", "OLS"))





