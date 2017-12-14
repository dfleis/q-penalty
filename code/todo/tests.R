mycols <- colorRampPalette(c("darkblue", "blue", "turquoise", "yellow", "orange")) # parula
#### functions 
blasso <- function(b, lam) { # orthogonal lasso function
  max(b - lam, 0)
}
bridge <- function(b, lam) { # orthogonal ridge function
  b/(1 + lam) # is it b/(2 * (1 + lam))
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
bqHALF_complex <- function(b_in, lam_in, upper = T, eps = 0) {
  # permit complex inputs for b_in
  bOLSmin <- (lam/4)^(2/3) + lam/2 * ((lam/4)^(2/3))^(-1/2)
  
  #bsign <- sign(b_in)
  lam <- complex(real = lam_in)
  I <- complex(real = 0, imaginary = 1)
  b <- b_in
  #b <- complex(real = b_in)
  
  term1 <-  2 * b / 3
  term2 <- -(1 - I * sqrt(3)) * b^2 / (3 * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3))
  term3 <- -1/12 * (1 + I * sqrt(3)) * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3)
  
  if (!upper) {
    term2 <- -(1 + I * sqrt(3)) * b^2 / (3 * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3))
    term3 <- -1/12 * (1 - I * sqrt(3)) * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3)
  }
  
  x <- term1 + term2 + term3
  
  return (x)
}


eps <- 0
lam <- 1
bmax_plot <- 3
bmin_plot <- 0

# administrative details
bmax <- bmax_plot * 1.1
bQmin <- (lam/4)^(2/3)
bOLSmin <- bQmin + lam/2 * (bQmin)^(-1/2)

# curves
beta_ols <- seq(0, bmax, length.out = 4e3)
beta_ols2 <- beta_ols # hacky solution to plot x ~ x
beta_q_up <- sapply(beta_ols, function(b) bqHALF(b_in = b, lam_in = lam, upper = T, eps = eps))
beta_q_lo <- sapply(beta_ols, function(b) bqHALF(b_in = b, lam_in = lam, upper = F, eps = eps))
beta_lasso <- sapply(beta_ols, function(b) blasso(b, lam))
beta_ridge <- sapply(beta_ols, function(b) bridge(b, lam))

# plots
plot(NA, 
     ylab = expression(beta),
     xlab = expression(beta^OLS),
     xlim = c(bmin_plot, bmax_plot), 
     ylim = c(0, max(Re(beta_q_up[!is.nan(beta_q_up)]))),
     main = substitute(
       paste("q = 0.5, ", lambda, " = ", lam, ", ", epsilon, " = ", eps), 
       list(lam = lam, eps = eps)))
abline(h = 0, col = 'gray80', lwd = 1, lty = 'dotted') # x axis
abline(v = 0, col = 'gray80', lwd = 1, lty = 'dotted') # y axis
lines(beta_ols ~ beta_ols2, lwd = 0.75, col = 'gray70')
lines(Re(beta_q_up) ~ beta_ols, lwd = 2, lty = 'solid', col = 'blue4')
lines(Re(beta_q_lo) ~ beta_ols, lwd = 2, lty = 'solid', col = 'deepskyblue')
points(bQmin ~ bOLSmin, pch = 19, cex = 0.25, col = 'blue4')
lines(beta_lasso ~ beta_ols, lwd = 1, col = 'darkgreen', lty = 'dashed')
lines(beta_ridge ~ beta_ols, lwd = 1, col = 'lightgreen', lty = 'dashed')
legend("topleft", col = c("blue4", "deepskyblue", "gray70", "darkgreen", "lightgreen"),
       lwd = c(2, 2, 0.75, 1, 1), bty = 'n', seg.len = 2, cex = 0.75,
       lty = c("solid", "solid", "solid", "dashed", "dashed"),
       legend = c("1/2-norm top branch", "1/2-norm top branch", "OLS", "Lasso", "Ridge"))



require(akima)
require(rgl)
lam <- 1
x <- seq(-2, 2, length.out = 500)
y <- seq(-2, 2, length.out = 500)
x_c <- complex(real = x)
y_c <- complex(imaginary = y)
Z <- outer(x_c, y_c, FUN = function(x, y) bqHALF_complex(b_in = x + y, lam_in = lam, upper = T, eps = 0))
z <- Re(Z)

surface3d(x, y, z, col = mycols(length(z)))
# pmat <- persp(x, y, z, phi = 30, theta = -50, 
#               xlim = c(min(x), max(x)),
#               ylim = c(min(y), max(y)),
#               xlab = "Re", ylab = "Im",
#               col = "springgreen", shade = 0.5)














