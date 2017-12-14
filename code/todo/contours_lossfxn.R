setwd("~/drive/projects/0_qnorm_penalized_regression/")

parula <- colorRampPalette(c(
  rgb(53, 42, 135, maxColorValue = 255),
  rgb(15, 92, 221, maxColorValue = 255),
  rgb(18, 125, 216, maxColorValue = 255),
  rgb(7, 156, 207, maxColorValue = 255),
  rgb(21, 177, 180, maxColorValue = 255),
  rgb(89, 189, 140, maxColorValue = 255),
  rgb(165, 190, 107, maxColorValue = 255),
  rgb(225, 185, 82, maxColorValue = 255),
  rgb(252, 206, 46, maxColorValue = 255),
  rgb(249, 251, 14, maxColorValue = 255)))

#### FUNCTIONS ####

L_fxn <- function(b, bols, lambda, q) {
  # orthogonal penalized loss function
  as.numeric(-t(b) %*% bols + 1/2 * t(b) %*% b + lambda * sum(abs(b)^q))  
}
L_fxn_deriv <- function(b, bols, lambda, q) {
  # orthogonal penalized loss function
  as.numeric(-bols + b + lambda * q * sum(abs(b)^(q - 1)))  
}

plotL_fxn <- function(b1, b2, bols, lambda, q, k, plotnorm = T) {
  bridge <- bols/(1 + 2 * lambda)
  blasso <- sign(bols) * pmax(abs(bols) - lambda, 0)
  
  L <- matrix(NA, nrow = length(b1), ncol = length(b2))
  
  for (i in 1:nrow(L)) {
    for (j in 1:ncol(L)) {
      L[i, j] <- L_fxn(c(b1[i], b2[j]), bols = bols, lambda = lambda, q = q)
    }
  }
  
  Lmin_idx <- which(L == min(L), arr.ind = TRUE)
  z <- contourLines(x = b1, y = b2, L, nlevels = 50)
  
  p <- image(b1, b2, L, col = NA, #parula(100),
        xlab = expression(beta[1]),
        ylab = expression(beta[2]),
        main = substitute(
          paste("q = ", q, ", ", lambda, " = ", lam), 
          list(q = q, lam = round(lambda, 3))))
  p <- p + abline(h = 0, col = "gray50", lty = 'dashed')
  p <- p + abline(v = 0, col = "gray50", lty = 'dashed')
  for (i in 1:length(z)) {
    p <- p + lines(z[[i]]$y ~ z[[i]]$x, lwd = 1, col = 'gray70')
  }
  
  if (plotnorm) {
    x <- seq(min(b1), max(b1), length.out = length(b1))
    y <- seq(min(b2), max(b2), length.out = length(b2))
  
    QN <- outer(x, y, FUN = function(x, y) (abs(x)^q + abs(y)^q)^(1/q))
    qnorm_lines <- contourLines(x = x, y = y, z = QN, nlevels = 1, levels = k)[[1]]
    p <- p + lines(qnorm_lines$y ~ qnorm_lines$x, lwd = 1.5, col = 'black', lty = 'solid')
  }
  
  p <- p + points(bols[1], bols[2], pch = 19, col = 'red', cex = 0.75)
  p <- p + text(bols[1], bols[2] * 1.1, expression(hat(beta)^OLS), cex = 0.75)
  
  if (q == 2) {
    p <- p + points(bridge[1], bridge[2], pch = 18, col = 'red', cex = 0.75)
    p <- p + text(bridge[1] * 1.4, bridge[2], expression(hat(beta)^ridge), cex = 0.75)
  } else if (q == 1) {
    p <- p + points(blasso[1], blasso[2], pch = 17, col = 'red', cex = 0.75)
    p <- p + text(blasso[1], blasso[2] * 0.9, expression(hat(beta)^lasso), cex = 0.75)
  } else if (lambda != 0) {
    p <- p + points(b1[Lmin_idx[1]], b2[Lmin_idx[2]], pch = 17, col = 'red', cex = 0.75)
    p <- p + text(b1[Lmin_idx[1]], b2[Lmin_idx[2]] * 0.9, expression(hat(beta)^q), cex = 0.75)
  }
}

#### 2-D Case ####
beta_ols <- c(0.8, 0.5)

b1min <- -0.1
b1max <- 1
b2min <- -0.1
b2max <- 1

b1 <- seq(b1min, b1max, length.out = 2e2)
b2 <- seq(b2min, b2max, length.out = 2e2)

lambdas <- seq(0, 1, length.out = 200)
q <- 0.5
k <- 0.5

plotL_fxn(b1, b2, beta_ols, lambda = lambdas[2], q = q, k = k, plotnorm = F)


max_z <- floor(log10(length(lambdas))) # max nb of leading zeros to add
pt <- proc.time()
for (idx in 1:length(lambdas)) {
  if (idx %% 5 == 0)
    print(idx)
  
  nb_digits <- floor(log10(idx)) + 1 # digits in our index
  nb_z <- max_z - nb_digits + 1 # nb of leading zeros to add
  idx_z <- paste0(c(rep(0, nb_z), idx), collapse = "")
  filename <- paste0(c("./img/plot/", "plot", idx_z, ".png"), collapse = "")
  
  png(filename, width = 600, height = 600)
  plotL_fxn(b1, b2, beta_ols, lambda = lambdas[idx], q = q, k = k, plotnorm = F)
  dev.off()
}
proc.time() - pt

plotdir <- paste0(c(getwd(), "/img/plot/"), collapse = "")
sys_cmd1 <- paste0("convert ", plotdir, "*.png -loop 0 ./img/gif/q", q, "_tmp.gif") # turn plots into tmp.gif
system(sys_cmd1)

sys_cmd2 <- paste0("convert -delay 12 ./img/gif/q", q, "_tmp.gif ./img/gif/q", q, "_out.gif") # increase tmp.gif framerate to 8ms/frame
system(sys_cmd2)

file.remove(list.files(plotdir, full.names = T)) # get full file patsh and remove them
list.files(plotdir)




#### 1-D Case ###

beta_ols <- 1.5
bmin <- -0.25
bmax <- 2
b <- seq(bmin, bmax, length.out = 1e3)

lambda <- 1
q <- 1

Lcurve <- sapply(b, L_fxn, bols = beta_ols, lambda = lambda, q = q)
Lcurve_deriv <- sapply(b, L_fxn_deriv, bols = beta_ols, lambda = lambda, q = q)
plot(NA,
     xlim = c(min(b), max(b)),
     ylim = c(min(Lcurve, Lcurve_deriv), 1),
     xlab = expression(beta),
     ylab = expression(L(beta)),
     main = substitute(
       paste("q = ", q, ", ", lambda, " = ", lam), 
       list(q = q, lam = lambda)))
abline(h = 0, col = "gray50", lty = 'dotted')
abline(v = 0, col = "gray50", lty = 'dotted')
lines(Lcurve ~ b, lwd = 2)
lines(Lcurve_deriv ~ b, lwd = 2, lty = 'dotted', col = 'red')



































