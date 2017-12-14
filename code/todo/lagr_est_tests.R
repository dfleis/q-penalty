##### setup #####
library(nloptr)
library(alabama)
library(NlcOptim)
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
mygrad <- colorRampPalette(c("red", "black"))

##### define some parameters #####
beta1 <- 0.8
beta2 <- 0.95
b1min <- 0
b1max <- 1
nb1 <- 100
b2min <- 0
b2max <- 1
nb2 <- 100

qs <- c(0.1 * (1:20), 3:10)
qs <- qs[which(qs != 0)]

lammin <- 0
lammax <- 20
nlams <- 500

b1 <- seq(b1min, b1max, length.out = nb1)
b2 <- seq(b2min, b2max, length.out = nb2)
lams <- seq(lammin, lammax, length.out = nlams)

##### define some functions #####
eps <- 0.001
obj_fxn <- function(x, y) -(x * beta1 + y * beta2) + 1/2 * (x^2 + y^2)
obj_fxn_lagr1 <- function(x, y, lam, q) obj_fxn(x, y) + lam * (abs(x)^q + abs(y)^q)
obj_fxn_lagr2 <- function(x, y, lam, q) obj_fxn(x, y) + lam * (sqrt(x^2 + eps^2)^q + sqrt(y^2 + eps^2)^q)

##### start work #####
loss_lagr_opt_list1 <- list()
loss_lagr_opt_list2 <- list()

# baseline loss contours/surface
loss_baseline <- outer(b1, b2, function(x, y) obj_fxn(x, y))

pt <- proc.time()
for (q_idx in 1:length(qs)) {
  q <- qs[q_idx]
  ## lagragian multiplier framework
  loss_lagr_opt_lams1 <- sapply(lams, FUN = function(l) 
    optim(c(0, 0), fn = function(params) obj_fxn_lagr1(params[1], params[2], l, q))$par)
  loss_lagr_opt_lams2 <- sapply(lams, FUN = function(l) 
    optim(c(0, 0), fn = function(params) obj_fxn_lagr2(params[1], params[2], l, q))$par)
  
  loss_lagr_opt_list1[[q_idx]] <- loss_lagr_opt_lams1
  loss_lagr_opt_list2[[q_idx]] <- loss_lagr_opt_lams2
}
proc.time() - pt

##### plots #####
## lagragian multiplier
# abs(x)
image(b1, b2, loss_baseline, col = parula(1e2),
      xlab = expression(beta[1]),
      ylab = expression(beta[2]),
      main = paste("lagr opt 1"))
for (i in 1:length(qs)) {
  mylwd <- 1
  mycol <- mygrad(length(qs))[i]
  if ((qs[i] == 1) || (qs[i] == 2)) {
    mycol <- "black"
    mylwd <- 2.5
  }
  lines(loss_lagr_opt_list1[[i]][2,] ~ loss_lagr_opt_list1[[i]][1,], col = mycol,
        type = 'l', cex = 0.75, pch = 19, lwd = mylwd)
}
# sqrt(x^2 + eps)
image(b1, b2, loss_baseline, col = parula(1e2),
      xlab = expression(beta[1]),
      ylab = expression(beta[2]),
      main = paste("lagr opt 2"))
#contour(b1, b2, loss_baseline, nlevels = 10, drawlabels = F, add = T)
for (i in 1:length(qs)) {
  mylwd <- 1
  mycol <- mygrad(length(qs))[i]
  if ((qs[i] == 1) || (qs[i] == 2)) {
    mycol <- "black"
    mylwd <- 2.5
  }
  lines(loss_lagr_opt_list2[[i]][2,] ~ loss_lagr_opt_list2[[i]][1,], col = mycol,
        type = 'l', cex = 0.75, pch = 19, lwd = mylwd)
}














