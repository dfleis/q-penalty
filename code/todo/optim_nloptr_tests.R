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

##### define some parameters #####
beta1 <- 0.7
beta2 <- 0.95
b1min <- 0
b1max <- 1
nb1 <- 150
b2min <- 0
b2max <- 1
nb2 <- 150

q <- 0.25

lammin <- 0
lammax <- 2
nlams <- 1e3
kmin <- 1
kmax <- 3
nks <- 1e3

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

##### start work #####
## baseline loss contours/surface
loss_baseline <- outer(b1, b2, function(x, y) obj_fxn(x, y))

## lagragian multiplier framework
pt <- proc.time()
loss_lagr_opt_lams <- sapply(lams, FUN = function(l) 
  optim(c(0, 0), fn = function(params) obj_fxn_lagr(params[1], params[2], l, q))$par)
proc.time() - pt

## constrained framework
# auglag: augmented lagragian algorithm
pt <- proc.time()
loss_constr_opt_ks <- sapply(ks, function(k)
  alabama::auglag(par = c(0, 0),
                  fn  = function(params) obj_fxn(params[1], params[2]),
                  gr  = function(params) obj_grad_fxn(params[1], params[2]),
                  hin = function(params) constr_fxn(params[1], params[2], k, q),
                  control.outer = list(maxeps = 1e-4, itmax = 1e4, trace = F))$par)
proc.time() - pt

# # COBYLA: constrained optim by linear approx
# pt <- proc.time()
# loss_constr_opt_ks2 <- sapply(ks, function(k)
#  cobyla(x0      = c(0, 0),
#         fn      = function(params) obj_fxn(params[1], params[2]),
#         lower   = c(0, 0),
#         upper   = c(beta1, beta2),
#         hin     = function(params) constr_fxn(params[1], params[2], k, q))$par)
# proc.time() - pt

# # ISRES: improved stochastic ranking evolution strategy
# pt <- proc.time()
# loss_constr_opt_ks3 <- sapply(ks, function(k)
#   nloptr(x0          = c(0, 0),
#          eval_f      = function(params) obj_fxn(params[1], params[2]),
#          lb          = c(0, 0),
#          ub          = c(beta1, beta2),
#          eval_g_ineq = function(params) constr_fxn2(params[1], params[2], k, q),
#          opts = list("algorithm" = "NLOPT_GN_ISRES",
#                      "xtol_rel"  = 1e-4,
#                      "maxeval"   = 1e4))$solution)
# proc.time() - pt

##### plots #####
# lagragian multiplier
image(b1, b2, loss_baseline, col = parula(1e2),
      xlab = expression(beta[1]),
      ylab = expression(beta[2]),
      main = substitute(paste("langr mult: q = ", q),
                        list(q = q)))
contour(b1, b2, loss_baseline, nlevels = 10, drawlabels = F, add = T)
lines(loss_lagr_opt_lams[2,] ~ loss_lagr_opt_lams[1,], col = 'red')
points(loss_lagr_opt_lams[2, 1] ~ loss_lagr_opt_lams[1, 1],
       pch = 19, cex = 0.75, col = 'red')

# auglag
image(b1, b2, loss_baseline, col = parula(1e2),
      xlab = expression(beta[1]),
      ylab = expression(beta[2]),
      main = substitute(paste("constr opt: auglag, q = ", q),
                        list(q = q)))
contour(b1, b2, loss_baseline, nlevels = 10, drawlabels = F, add = T)
lines(loss_constr_opt_ks[2,] ~ loss_constr_opt_ks[1,], col = 'red')
points(rev(loss_constr_opt_ks[2, ])[1] ~ rev(loss_constr_opt_ks[1, ])[1], 
       pch = 19, cex = 0.5, col = 'red')

# # cobyla
# image(b1, b2, loss_baseline, col = parula(1e2),
#       xlab = expression(beta[1]),
#       ylab = expression(beta[2]),
#       main = substitute(paste("constr opt: cobyla, q = ", q),
#                         list(q = q)))
# contour(b1, b2, loss_baseline, nlevels = 10, drawlabels = F, add = T)
# lines(loss_constr_opt_ks2[2,] ~ loss_constr_opt_ks2[1,], col = 'red')
# points(rev(loss_constr_opt_ks2[2, ])[1] ~ rev(loss_constr_opt_ks2[1, ])[1],
#        pch = 19, cex = 0.75, col = 'red')

# # isres
# image(b1, b2, loss_baseline, col = parula(1e2),
#       xlab = expression(beta[1]),
#       ylab = expression(beta[2]),
#       main = substitute(paste("constr opt: isres q = ", q),
#                         list(q = q)))
# contour(b1, b2, loss_baseline, nlevels = 10, drawlabels = F, add = T)
# lines(loss_constr_opt_ks3[2,] ~ loss_constr_opt_ks3[1,], col = 'red')
# points(rev(loss_constr_opt_ks3[2, ])[1] ~ rev(loss_constr_opt_ks3[1, ])[1],
#        pch = 19, cex = 0.75, col = 'red')









