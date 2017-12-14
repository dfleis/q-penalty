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
betas <- c(0.7, 0.95, 0.88)
bmins <- c(0, 0, 0)
bmaxs <- c(1, 1, 1)
nbs <- c(10, 10, 10)

qs <- 0.5
qs <- qs[which(qs != 0)]

kmin <- 0
kmax <- 3
nks <- 10

bparams <- cbind(bmins, bmaxs, nbs)
bparams <- split(bparams, seq(ncol(bparams))) # turn matrix into a list by rows of the matrix
bs <- lapply(bparams, FUN = function(bpar) seq(from = bpar[1], to = bpar[2], length.out = bpar[3]))

lams <- 1#seq(lammin, lammax, length.out = nlams)
ks <- seq(kmin, kmax, length.out = nks)

##### define some functions #####
obj_fxn <- function(v) -as.numeric(v %*% betas + 1/2 * v %*% v)
obj_grad_fxn <- function(v) -betas + v
constr_fxn <- function(v, k, q) as.numeric(k - sum(abs(v)^q)) # hin(x) >= 0
constr_fxn2 <- function(v, k, q) as.numeric(sum(abs(v)^q) - k) # g_ineq(x) <= 0
obj_fxn_lagr <- function(v, lam, q) as.numeric(obj_fxn(v) + lam * sum(abs(v)^q))

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













