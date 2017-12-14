mycols <- colorRampPalette(c("darkblue", "blue", "turquoise", "yellow", "orange")) # parula

no_soln_fxn <- function(q, lam) {
  (lam * (1 - q) * q)^(1/(2 - q)) + lam * q * ( (lam * (1 - q) * q)^(1/(2 - q)) )^(q - 1)
}

mylams <- c(1, 3/4, 1/2, 1/4)
cols <- c("black", "red", "blue", "green3")
ltys <- c("solid", "dashed", "dotted", "dotdash")
q <- seq(0.001, 1, by = 0.001)

plot(no_soln_fxn(q, mylams[1]) ~ q, type = 'l',
     col = cols[1], lty = ltys[1],
     xlab = "q", ylab = "Width",
     main = "Interval Width")
for (i in 2:length(mylams))
  lines(no_soln_fxn(q, mylams[i]) ~ q, col = cols[i], lty = ltys[i])
legend("topleft", legend = 
         c(eval(substitute(expression(paste(lambda, " = ", l)), list(l = mylams[1]))), 
           eval(substitute(expression(paste(lambda, " = ", l)), list(l = mylams[2]))), 
           eval(substitute(expression(paste(lambda, " = ", l)), list(l = mylams[3]))),
           eval(substitute(expression(paste(lambda, " = ", l)), list(l = mylams[4])))),
       col = cols, seg.len = 2, cex = 0.75, lwd = 2,
       lty = ltys, bty = 'n')







q <- seq(0.001, 1, by = 0.0025)
lam <- seq(0.1, 10, length.out = 2.5e2)
max_no_soln <- outer(q, lam, FUN = function(q, lam) no_soln_fxn(q, lam))

contour(q, lam, max_no_soln, xlab = "q", ylab = expression(lambda))
image(q, lam, max_no_soln, xlab = "q", ylab = expression(lambda), col = mycols(100))

library(fields)
image.plot(q, lam, max_no_soln, xlab = "q", ylab = expression(lambda), col = mycols(100))
contour(q, lam, max_no_soln, xlab = "q", ylab = expression(lambda), add = T, lwd = 2)
