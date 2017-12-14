


bPOS <- function(x, q, lam) {
  x + lam * q * x^(q - 1)
}
bNEG <- function(x, q, lam) {
  x - lam * q * (-x)^(q - 1)
}
bPOS2 <- function(x, q, lam) {
  betaQmin <- (q * lam - q^2 * lam)^(1/(2 - q))
  
  if (x < betaQmin) {
    return (NaN)
  } else {
    x + lam * q * x^(q - 1)
  }
}
bNEG2 <- function(x, q, lam) {
  betaQmin <- (q * lam - q^2 * lam)^(1/(2 - q))
  
  if (abs(x) < betaQmin) {
    return (NaN)
  } else {
    x - lam * q * (-x)^(q - 1)
  }
}

blasso <- function(bols, lam) {
  sign(bols) * max(abs(bols) - lam, 0)
}
bqHALF <- function(b, lam) {
  bsign <- sign(b)
  b <- complex(real = abs(b))
  lam <- complex(real = lam)
  I <- complex(real = 0, imaginary = 1)
  
  term1 <-  2 * b / 3
  term2 <- -(1 - I * sqrt(3)) * b^2 / (3 * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3))
  term3 <- -1/12 * (1 + I * sqrt(3)) * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3)
  
  x <- term1 + term2 + term3
  
  if (abs(Im(x)) > 1e-3) {
    return (0)
  } else {
    return (bsign * x)
  }
}


q <- 1/2
lam <- 4/ (3 * sqrt(3))
lam <- 1
bmax <- 4


#### PLOT ALL BRANCHES
betaols <- seq(0.01, bmax, length.out = 1e2)
bP <- sapply(betaols, function(z) bPOS(z, q, lam))
bN <- sapply(-betaols, function(z) bNEG(z, q, lam))

betaols_input <- seq(-max(bP, na.rm = T), max(bP, na.rm = T), length.out = 1e2)
bL <- sapply(betaols_input, function(z) blasso(z, lam))
plot(betaols ~ bP, type = 'l', lwd = 2, 
     xlim = c(-max(bP, na.rm = T), max(bP, na.rm = T)),
     ylim = c(-max(bP, na.rm = T), max(bP, na.rm = T)),
     xlab = expression(beta^OLS),
     ylab = expression(beta),
     main = substitute(
       paste("q = ", q, ", ", lambda, " = ", lam),
       list(q = q, lam = lam)))
lines(-betaols ~ bN, lwd = 2)
abline(h = 0, col = 'gray50', lty = 'dashed')
abline(v = 0, col = 'gray50', lty = 'dashed')
lines(bL ~ betaols_input, lwd = 1.5, col = "blue")

bQH <- sapply(betaols_input, function(b) bqHALF(b, lam))
lines(Re(bQH) ~ betaols_input, type = 'l', col = 'red', lwd = 2, lty = 'dashed')

betaQmin <- (q * lam - q^2 * lam)^(1/(2 - q))



betaOLSmin <- betaQmin + lam * q * betaQmin^(q - 1)
abline(h = betaQmin, col = 'green3', lty = 'dotted')
abline(v = betaOLSmin, col = 'green3', lty = 'dotted')
abline(h = -betaQmin, col = 'green3', lty = 'dotted')
abline(v = -betaOLSmin, col = 'green3', lty = 'dotted')

####### only plot extreme branches
betaols <- seq(0.01, bmax, length.out = 1e2)
bP2 <- sapply(betaols, function(z) bPOS2(z, q, lam))
bN2 <- sapply(-betaols, function(z) bNEG2(z, q, lam))

betaols_input <- seq(-max(bP, na.rm = T), max(bP, na.rm = T), length.out = 1e2)
bL <- sapply(betaols_input, function(z) blasso(z, lam))
plot(betaols ~ bP2, type = 'l', lwd = 2, 
     xlim = c(-max(bP, na.rm = T), max(bP, na.rm = T)),
     ylim = c(-max(bP, na.rm = T), max(bP, na.rm = T)),
     xlab = expression(beta^OLS),
     ylab = expression(beta),
     main = substitute(
       paste("q = ", q, ", ", lambda, " = ", lam),
       list(q = q, lam = lam)))
lines(-betaols ~ bN2, lwd = 2)
abline(h = 0, col = 'gray50', lty = 'dashed')
abline(v = 0, col = 'gray50', lty = 'dashed')
lines(bL ~ betaols_input, lwd = 1.5, col = "blue")

bQH <- sapply(betaols_input, function(b) bqHALF(b, lam))
lines(Re(bQH) ~ betaols_input, type = 'l', col = 'red', lwd = 2, lty = 'dashed')




