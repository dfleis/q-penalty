blasso <- function(bols, lam) {
  sign(bols) * max(abs(bols) - lam, 0)
}
bQTop <- function(b, lam) {
  bsign <- sign(Re(b))
  b <- complex(real = abs(b))
  lam <- complex(real = lam)
  I <- complex(real = 0, imaginary = 1)
  
  term1 <-  2 * b / 3
  term2 <- -(1 - I * sqrt(3)) * b^2 / (3 * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3))
  term3 <- -1/12 * (1 + I * sqrt(3)) * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3)
  
  x <- term1 + term2 + term3
  
  return (bsign * x)
}
bQBottom <- function(b, lam) {
  bsign <- sign(Re(b))
  b <- complex(real = abs(b))
  lam <- complex(real = lam)
  I <- complex(real = 0, imaginary = 1)
  
  term1 <-  2 * b / 3
  term2 <- -(1 + I * sqrt(3)) * b^2 / (3 * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3))
  term3 <- -1/12 * (1 - I * sqrt(3)) * (-8 * b^3 + 27 * lam^2 + 3 * sqrt(3) * sqrt(-16 * b^3 * lam^2 + 27 * lam^4))^(1/3)
  
  x <- term1 + term2 + term3

  return (bsign * x)
}




lam <- 1
bmax <- 2
b <- seq(-bmax, bmax, length.out = 1e3)
b2 <- b # hacky solution to plot x ~ x

bQT <- sapply(b, FUN = function(b) bQTop(b, lam))
bQB <- sapply(b, FUN = function(b) bQBottom(b, lam))
bL <- sapply(b, FUN = function(b) blasso(b, lam))

plot(abs(bQT) ~ b, type = 'l',
     xlab = expression(beta^OLS),
     ylab = expression(beta),
     xlim = c(-bmax, bmax),
     ylim = c(-max(Re(bQT), na.rm = T), max(Re(bQT), na.rm = T)),
     lwd = 2, col = 'blue',
     main = substitute(
       paste(lambda, " = ", lam),
       list(lam = lam)))
lines(Re(bQB) ~ b, lwd = 2, col = 'red')
abline(h = 0, col = 'gray50', lty = 'dashed')
abline(v = 0, col = 'gray50', lty = 'dashed')



plot(bQT, pch = 19, cex = 0.25)
plot(bQB, pch = 19, cex = 0.25)

re <- seq(0, 2, length.out = 1e2)
im <- complex(real = rep(0, length = length(im)), imaginary = re)
M <- outer(re, im, FUN = function(a, b) {
  z <- a + b
  bQTop(z, lam)
})

plot(abs(M))






