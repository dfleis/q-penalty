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
beta1 <- 0.7
beta2 <- 0.95
b1min <- 0
b1max <- 1
nb1 <- 100
b2min <- 0
b2max <- 1
nb2 <- 100

qs <- 0.1 * (1:10)
qs <- qs[which(qs != 0)]

lammin <- 0
lammax <- 10
nlams <- 1e2
kmin <- 0
kmax <- 3
nks <- 1e2

b1 <- seq(b1min, b1max, length.out = nb1)
b2 <- seq(b2min, b2max, length.out = nb2)
lams <- seq(lammin, lammax, length.out = nlams)
ks <- seq(kmin, kmax, length.out = nks)

##### define some functions #####
obj_fxn <- function(x, y) -(x * beta1 + y * beta2) + 1/2 * (x^2 + y^2)
obj_grad_fxn <- function(x, y) c(-beta1 + x, -beta2 + y)
constr_fxn <- function(x, y, k, q) k - abs(x)^q - abs(y)^q # hin(x) >= 0
constr_fxn2 <- function(x, y, k, q) abs(x)^q + abs(y)^q - k # g_ineq(x) <= 0
obj_fxn_lagr <- function(x, y, lam, q) obj_fxn(x, y) + lam * (abs(x)^q + abs(y)^q)
eps <- 0.1
obj_fxn_lagr <- function(x, y, lam, q) obj_fxn(x, y) + lam * (sqrt(x^2 + eps)^q + sqrt(y^2 + eps)^q)

##### start work #####
loss_constr_opt_list <- list()
loss_lagr_opt_list <- list()

# baseline loss contours/surface
loss_baseline <- outer(b1, b2, function(x, y) obj_fxn(x, y))

pt <- proc.time()
for (q_idx in 1:length(qs)) {
  q <- qs[q_idx]
  ## lagragian multiplier framework
  loss_lagr_opt_lams <- sapply(lams, FUN = function(l) 
    optim(c(0, 0), fn = function(params) obj_fxn_lagr(params[1], params[2], l, q))$par)
  
  ## constrained framework
  # auglag: augmented lagragian algorithm
  loss_constr_opt_ks <- sapply(ks, function(k)
    alabama::auglag(par = c(0, 0),
                    fn  = function(params) obj_fxn(params[1], params[2]),
                    gr  = function(params) obj_grad_fxn(params[1], params[2]),
                    hin = function(params) constr_fxn(params[1], params[2], k, q),
                    control.outer = list(maxeps = 1e-4, itmax = 1e4, trace = F))$par)
  
  loss_lagr_opt_list[[q_idx]] <- loss_lagr_opt_lams
  loss_constr_opt_list[[q_idx]] <- loss_constr_opt_ks
}
proc.time() - pt

##### plots #####
# lagragian multiplier
image(b1, b2, loss_baseline, col = parula(1e2),
      xlab = expression(beta[1]),
      ylab = expression(beta[2]),
      main = paste("lagr opt"))
#contour(b1, b2, loss_baseline, nlevels = 10, drawlabels = F, add = T)
for (i in 1:length(qs)) {
  mylwd <- 1
  mycol <- mygrad(length(qs))[i]
  if ((qs[i] == 1) || (qs[i] == 2)) {
    mycol <- "black"
    mylwd <- 2.5
  }
  lines(loss_lagr_opt_list[[i]][2,] ~ loss_lagr_opt_list[[i]][1,], col = mycol,
        type = 'l', cex = 0.75, pch = 19, lwd = mylwd)
}

# auglag
image(b1, b2, loss_baseline, col = parula(1e2),
      xlab = expression(beta[1]),
      ylab = expression(beta[2]),
      main = paste("constr opt: auglag"))
#contour(b1, b2, loss_baseline, nlevels = 10, drawlabels = F, add = T)
for (i in 1:length(qs)) {
  mylwd <- 1
  mycol <- mygrad(length(qs))[i]
  if ((qs[i] == 1) || (qs[i] == 2)) {
    mycol <- "black"
    mylwd <- 2.5
  }
  lines(loss_constr_opt_list[[i]][2,] ~ loss_constr_opt_list[[i]][1,], col = mycol,
        type = 'l', cex = 0.75, pch = 19, lwd = mylwd)
}













