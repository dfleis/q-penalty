mycols <- colorRampPalette(c("black", "green2"))


bPOS <- function(x, q, lam) {
  x + lam * q * x^(q - 1)
}
bPOS2 <- function(x, q, lam) {
  betaQmin <- (q * lam - q^2 * lam)^(1/(2 - q))
  
  if (x < betaQmin) {
    return (NaN)
  } else {
    x + lam * q * x^(q - 1)
  }
}
bPOS3 <- function(x, q, lam) {
  betaQmin <- (q * lam - q^2 * lam)^(1/(2 - q))
  
  if (x > betaQmin) {
    return (NaN)
  } else {
    x + lam * q * x^(q - 1)
  }
}

blasso <- function(bols, lam) {
  sign(bols) * max(abs(bols) - lam, 0)
}


q <- 1/2
lam <- 1
bmax <- 4


#### PLOT ALL BRANCHES
betaols <- seq(0.01, bmax, length.out = 1e3)
bP <- sapply(betaols, function(z) bPOS(z, q, lam))

betaols_input <- seq(0, max(bP, na.rm = T), length.out = 1e3)
betaols_input2 <- betaols_input # hacky solution to plot x ~ x
bL <- sapply(betaols_input, function(z) blasso(z, lam))
plot(betaols ~ bP, type = 'l', lwd = 2, 
     xlim = c(0, min(bmax, max(bP, na.rm = T))),
     ylim = c(0, min(bmax, max(bP, na.rm = T))),
     xlab = expression(beta^OLS),
     ylab = expression(beta),
     main = substitute(
       paste("q = ", q, ", ", lambda, " = ", lam),
       list(q = q, lam = lam)))
abline(h = 0, col = 'gray80')
abline(v = 0, col = 'gray80')
lines(bL ~ betaols_input, lwd = 1.5, col = 'gray50')
lines(betaols_input2 ~ betaols_input, lwd = 1.5, col = 'gray50', lty = 'dotted')
legend("topleft", legend = c("q Penalty", "Lasso", "OLS"), bty = 'n',
       lty = c("solid", "dashed", "dotted"), lwd = c(2, 1.5, 1.5), seg.len = 2.5,
       col = c("black", "gray50", "gray50"))




#### PLOT ALL BRANCHES as a function of lambda
q <- 1/2
lams <- c(1, 3/4, 2/4, 1/4)
bmax <- 4

cols <- c("black", "red", "blue", "green3")
ltys <- c("solid", "dashed", "dotted", "dotdash")


betaols <- seq(0.01, bmax, length.out = 1e3)
bP <- lapply(lams, FUN = function(lam) sapply(betaols, function(z) bPOS(z, q, lam)))

betaols_input <- seq(0, max(unlist(bP), na.rm = T), length.out = 1e2)
betaols_input2 <- betaols_input # hacky solution to plot x ~ x
bL <- lapply(lams, FUN = function(lam) sapply(betaols_input, function(z) blasso(z, lam)))
plot(betaols ~ bP[[1]], type = 'l', lwd = 2, 
     col = cols[[1]], lty = ltys[[1]],
     xlim = c(0, min(unlist(lapply(bP, max)))),
     ylim = c(0, min(unlist(lapply(bP, max)))),
     xlab = expression(beta^OLS),
     ylab = expression(beta),
     main = substitute(
       paste("q = ", q),
       list(q = q)))
for (i in 2:length(lams)) {
  lines(betaols ~ bP[[i]], lwd = 2, col = cols[i], lty = ltys[i])
}
abline(h = 0, col = 'gray80')
abline(v = 0, col = 'gray80')
lines(betaols_input2 ~ betaols_input, lwd = 1, col = 'gray80')


#### PLOT ALL BRANCHES as a function of q
lam <- 1
qs <- c(0.95, seq(0.8, 0.2, by = -0.2), 0.01)
bmax <- 4

cols <- mycols(length(qs))
ltys <- c("solid", "dashed", "dotted", "dotdash")
ltys <- rep("solid", length(qs))


betaols <- seq(0.000000000000001, bmax, length.out = 1e3)
bP <- lapply(qs, FUN = function(q) sapply(betaols, function(z) bPOS(z, q, lam)))

betaols_input <- seq(0, max(unlist(bP), na.rm = T), length.out = 1e2)
betaols_input2 <- betaols_input # hacky solution to plot x ~ x
bL <- sapply(betaols, function(z) blasso(z, lam))

plot(NA,
     xlim = c(0, min(unlist(lapply(bP, max)), bmax)),
     ylim = c(0, min(unlist(lapply(bP, max)), bmax)),
     xlab = expression(beta^OLS),
     ylab = expression(beta),
     main = substitute(
       paste(lambda, " = ", lam),
       list(lam = lam)))
abline(h = 0, col = 'gray80')
abline(v = 0, col = 'gray80')
lines(betaols_input2 ~ betaols_input, lwd = 1.5, col = 'gray50', lty = 'dotted')
for (i in (length(qs) - 1):2) {
  lines(betaols ~ bP[[i]], lwd = 1, col = cols[i], lty = ltys[i])
}
lines(betaols ~ bP[[length(qs)]], lwd = 2, col = cols[length(qs)], lty = ltys[length(qs)])
lines(bL ~ betaols, lwd = 1.5, col = 'red')
lines(betaols ~ bP[[1]], type = 'l', lwd = 2, col = cols[1], lty = ltys[1])

legend("topleft", legend = 
         c(eval(substitute(expression(paste("q = ", q)), list(q = qs[1]))), 
           eval(substitute(expression(paste("q = ", q)), list(q = qs[length(qs)]))),
           "OLS", "Lasso"),
       col = c(cols[1], cols[length(qs)], "red", "gray50"), seg.len = 2.5, cex = 0.75, lwd = 2,
       lty = c(ltys[1], ltys[length(qs)], "solid", "dotted"), bty = 'n')





