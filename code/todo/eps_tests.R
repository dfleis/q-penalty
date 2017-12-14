plotmax <- 2
q <- 1/2
lam <- 1
eps <- 0.1

x <- seq(-0.25, 2, length.out = 1e3)
x_im <- complex(real = x, imaginary = eps)

x2 <- x # hacky solution to plotting x ~ x
b <- x + lam * q * x^(q-1)
b_im <- x_im + lam * q * x_im^(q - 1)

# plots
plot(b ~ x, type = 'l', 
     lwd = 2,
     ylim = c(-0.25, plotmax),
     xlim = c(-0.25, plotmax),
     main = substitute(
       paste("q = ", q, ", ", lambda, " = ", lam, ", ", epsilon, " = ", eps), 
       list(lam = lam, q = q, eps = eps)
       ))
lines(Re(b_im) ~ x, col = 'red', lwd = 1)
lines(abs(Im(b_im)) ~ x, col = 'red', lwd = 1, lty = 'dotted')
lines(x ~ x2, col = 'gray70', lty = 'dashed')
abline(h = 0, col = 'gray70', lty = 'dotted') # x axis
abline(v = 0, col = 'gray70', lty = 'dotted') # y axis

plot(x ~ b, type = 'l', 
     lwd = 2,
     ylim = c(-0.25, plotmax),
     xlim = c(-0.25, plotmax),
     main = substitute(
       paste("q = ", q, ", ", lambda, " = ", lam, ", ", epsilon, " = ", eps), 
       list(lam = lam, q = q, eps = eps)
     ))
lines(x ~ Re(b_im), col = 'red', lwd = 1)
lines(x ~ abs(Im(b_im)), col = 'red', lwd = 1, lty = 'dotted')
lines(x ~ x2, col = 'gray70', lty = 'dashed')
abline(h = 0, col = 'gray70', lty = 'dotted') # x axis
abline(v = 0, col = 'gray70', lty = 'dotted') # y axis




