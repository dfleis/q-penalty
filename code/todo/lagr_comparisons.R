setwd("~/drive/projects/0_qnorm_penalized_regression/")
parula <- colorRampPalette(c(
  rgb(53,  42,  135, maxColorValue = 255),
  rgb(15,  92,  221, maxColorValue = 255),
  rgb(18,  125, 216, maxColorValue = 255),
  rgb(7,   156, 207, maxColorValue = 255),
  rgb(21,  177, 180, maxColorValue = 255),
  rgb(89,  189, 140, maxColorValue = 255),
  rgb(165, 190, 107, maxColorValue = 255),
  rgb(225, 185, 82,  maxColorValue = 255),
  rgb(252, 206, 46,  maxColorValue = 255),
  rgb(249, 251, 14,  maxColorValue = 255)))

## define some parameters
q <- 0.1
beta_ols <- c(0.7, 0.9)

kmin <- 0
kmax <- 3.5
dk   <- 0.01
lammin <- 0
lammax <- 0.6
dlam   <- 0.0025

b1min <- 0
b1max <- 1
db1   <- 0.001
b2min <- 0
b2max <- 1
db2   <- 0.001

b1   <- seq(b1min, b1max, db1)
b2   <- seq(b2min, b2max, db2)
ks   <- seq(kmin, kmax, dk)
lams <- seq(lammin, lammax, dlam)

## define some functions
S <- function(x, y) -(x * beta_ols[1] + y * beta_ols[2]) + 1/2 * (x^2 + y^2)
SLagr <- function(x, y, lam, q) S(x, y) + lam * (abs(x)^q + abs(y)^q)
norm_q <- function(x, y, q) outer(x, y, function(x, y) (abs(x)^q + abs(y)^q)^(1/q))

## do work
L <- outer(b1, b2, FUN = function(x, y) S(x, y))

x <- vector(mode = 'numeric', length = length(ks))
y <- vector(mode = 'numeric', length = length(ks))
xLagr <- vector(mode = 'numeric', length = length(lams))
yLagr <- vector(mode = 'numeric', length = length(lams))
LLagr_list <- list()

pt <- proc.time()
for (i in 1:length(ks)) {
  if (i %% 25 == 0)
    print(paste0("k index % = ", round(i/length(ks)*100, 3)))
  
  z <- outer(b1, b2, FUN = function(x, y) (abs(x)^q + abs(y)^q)^(1/q))
  pos <- which(z <= ks[i], arr.ind = TRUE)
  minpos <- which(L[pos] == min(L[pos]))[1]
  minposx <- pos[minpos, 1]
  minposy <- pos[minpos, 2]

  x[i] <- b1[minposx]
  y[i] <- b2[minposy]
}
for (i in 1:length(lams)) {
  if (i %% 25 == 0)
    print(paste0("lam index % = ", round(i/length(ks)*100, 3)))
  
  LLagr_list[[i]] <- outer(b1, b2, FUN = function(x, y) SLagr(x, y, lams[i], q))
  minposLagr <- which(LLagr_list[[i]] == min(LLagr_list[[i]]), arr.ind = TRUE)
  
  xLagr[i] <- b1[minposLagr[1]]
  yLagr[i] <- b2[minposLagr[2]]
}
proc.time() - pt

plot(x ~ ks, type = 'l', col = 'red',
     ylim = c(min(x, y), max(x, y)), 
     xlab = "t",
     ylab = "Coefficient")
lines(y ~ ks, col = 'blue')

plot(xLagr ~ lams, type = 'l', col = 'red',
     ylim = c(min(xLagr, yLagr), max(xLagr, yLagr)), 
     xlab = expression(lambda),
     ylab = "Coefficient")
lines(yLagr ~ lams, col = 'blue')

plot(yLagr ~ xLagr, type = 'l', col = "blue", lwd = 2.5,
     xlim = c(min(b1), max(b1)),
     ylim = c(min(b2), max(b2)),
     xlab = expression(beta[1]),
     ylab = expression(beta[2]))
lines(y ~ x, col = "red", lwd = 1.25)
points(beta_ols[1], beta_ols[2], pch = 19)
#legend("bottomright", legend = c("Base Formulation", "Lagragian"), bty = 'n',
#       lwd = c(1.25, 2.5), col = c("red", "blue"), seg.len = 2.5,
#       lty = c("solid", "solid"))

# gif 1
pt <- proc.time()
max_z <- ceiling(log10(length(lams)))
for (idx in 1:length(lams)) {
  if (idx %% 5 == 0)
    print(idx)
  
  nb_digits <- floor(log10(idx)) + 1 # digits in our index
  nb_z <- max_z - nb_digits + 1 # nb of leading zeros to add
  idx_z <- paste0(c(rep(0, nb_z), idx), collapse = "")
  filename <- paste0(c("./img/tmp_plot/", "plot", idx_z, ".png"), collapse = "")
  
  png(filename, width = 600, height = 600)
  image(LLagr_list[[idx]], col = parula(100),
        xlab = expression(beta[1]),
        ylab = expression(beta[2]),
        main = substitute(paste(lambda, " = ", l, ", q = ", q), 
                          list(l = formatC(lams[idx], format = "f", flag = '0', digits = 3),
                               q = q)))
  lines(yLagr[1:idx] ~ xLagr[1:idx], lwd = 2)
  points(yLagr[idx] ~ xLagr[idx], pch = 19, cex = 2)
  dev.off()
}
proc.time() - pt

plotdir <- paste0(c(getwd(), "/img/tmp_plot/"), collapse = "")
sys_cmd1 <- paste0("convert ", plotdir, "*.png -loop 0 ./img/gif/lagr_q", q, "_tmp.gif") # turn plots into tmp.gif
system(sys_cmd1)

sys_cmd2 <- paste0("convert -delay 10 ./img/gif/lagr_q", q, "_tmp.gif ./img/gif/lagr_q", q, "_out.gif") # increase tmp.gif framerate to 8ms/frame
system(sys_cmd2)

file.remove(list.files(plotdir, full.names = T)) # get full file patsh and remove them
list.files(plotdir)



## gif 2
pt <- proc.time()
max_z <- ceiling(log10(length(ks)))
for (idx in 1:length(ks)) {
  if (idx %% 5 == 0)
    print(idx)
  
  nb_digits <- floor(log10(idx)) + 1 # digits in our index
  nb_z <- max_z - nb_digits + 1 # nb of leading zeros to add
  idx_z <- paste0(c(rep(0, nb_z), idx), collapse = "")
  filename <- paste0(c("./img/tmp_plot/", "plot", idx_z, ".png"), collapse = "")
  
  png(filename, width = 600, height = 600)
  image(L, col = parula(100),
        xlab = expression(beta[1]),
        ylab = expression(beta[2]),
        main = substitute(paste(k, " = ", kk, ", q = ", q), 
                          list(kk = formatC(rev(ks)[idx], format = "f", flag = '0', digits = 3),
                               q = q)))
  contour(b1, b2, norm_q(b1, b2, q), nlevels = 1, levels = rev(ks)[idx], add = T,
          bty = 'n', xaxt = 'n', yaxt = 'n', lwd = 1.5, col = 'red')
  lines(rev(y)[1:idx] ~ rev(x)[1:idx], lwd = 2)
  points(rev(y)[idx] ~ rev(x)[idx], pch = 19, cex = 2)
  dev.off()
}
proc.time() - pt

plotdir <- paste0(c(getwd(), "/img/tmp_plot/"), collapse = "")
sys_cmd1 <- paste0("convert ", plotdir, "*.png -loop 0 ./img/gif/base_q", q, "_tmp.gif") # turn plots into tmp.gif
system(sys_cmd1)

sys_cmd2 <- paste0("convert -delay 10 ./img/gif/base_q", q, "_tmp.gif ./img/gif/base_q", q, "_out.gif") # increase tmp.gif framerate to 8ms/frame
system(sys_cmd2)

file.remove(list.files(plotdir, full.names = T)) # get full file patsh and remove them
list.files(plotdir)



