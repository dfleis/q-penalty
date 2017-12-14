


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


q <- 1/2
lam <- 4/ (3 * sqrt(3))
lam <- 1
bmax <- 4


#### PLOT ALL BRANCHES
betaols <- seq(0.01, bmax, length.out = 1e2)
bP <- sapply(betaols, function(z) bPOS(z, q, lam))
bN <- sapply(-betaols, function(z) bNEG(z, q, lam))

betaols_input <- seq(0, max(bP, na.rm = T), length.out = 1e2)
betaols_input2 <- betaols_input
bL <- sapply(betaols_input, function(z) blasso(z, lam))
plot(betaols ~ bP, type = 'l', lwd = 2, 
     xlim = c(0, max(bP, na.rm = T)),
     ylim = c(0, max(bP, na.rm = T)),
     xlab = expression(beta^OLS),
     ylab = expression(beta),
     main = substitute(
       paste("q = ", q, ", ", lambda, " = ", lam),
       list(q = q, lam = lam)))
lines(-betaols ~ bN, lwd = 2)
abline(h = 0, col = 'gray80', lty = 'dashed')
abline(v = 0, col = 'gray80', lty = 'dashed')
lines(bL ~ betaols_input, lwd = 1.5, col = 'gray50')
lines(betaols_input2 ~ betaols_input, lwd = 1.5, col = 'gray50')

