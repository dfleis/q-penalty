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
beta <- 0.7
bmin <- 0
bmax <- 1
nb <- 200

qs <- unique(c(10^-(3:2), seq(0, 1, length.out = 10), 
               seq(1, 2, length.out = 10),
               seq(2, 10, length.out = 10)))
qs <- qs[which(qs != 0)]

lammin <- 0
lammax <- 2
nlams <- 200
kmin <- 0
kmax <- 1
nks <- 200

b <- seq(bmin, bmax, length.out = nb)
lams <- seq(lammin, lammax, length.out = nlams)
ks <- seq(kmin, kmax, length.out = nks)

##### define some functions #####
obj_fxn <- function(x) -x * beta1 + 1/2 * x^2
obj_grad_fxn <- function(x) -beta + x
constr_fxn <- function(x, k, q) k - abs(x)^q # hin(x) >= 0
constr_fxn2 <- function(x, k, q) abs(x)^q - k # g_ineq(x) <= 0
obj_fxn_lagr <- function(x, lam, q) obj_fxn(x) + lam * abs(x)^q

##### start work #####
loss_constr_opt_list <- list()
loss_lagr_opt_list <- list()

# baseline loss contours/surface
loss_baseline <- sapply(b, obj_fxn)

pt <- proc.time()
for (q_idx in 1:length(qs)) {
  q <- qs[q_idx]
  ## lagragian multiplier framework
  loss_lagr_opt_lams <- sapply(lams, FUN = function(l) 
    optimize(f = function(params) obj_fxn_lagr(params, l, q), lower = 0, upper = beta1)$minimum)
  
  ## constrained framework
  # auglag: augmented lagragian algorithm
  loss_constr_opt_ks <- sapply(ks, function(k)
    alabama::auglag(par = c(0),
                    fn  = function(params) obj_fxn(params),
                    gr  = function(params) obj_grad_fxn(params),
                    hin = function(params) constr_fxn(params, k, q),
                    control.outer = list(maxeps = 1e-4, itmax = 1e4, trace = F))$par)
  
  loss_lagr_opt_list[[q_idx]] <- loss_lagr_opt_lams
  loss_constr_opt_list[[q_idx]] <- loss_constr_opt_ks
}
proc.time() - pt

##### plots #####
# lagragian multiplier
plot(loss_lagr_opt_list[[1]] ~ rev(lams), pch = 19, cex = 0.5, col = mygrad(length(qs))[1], type = 'l')
for (i in 2:length(qs)) {
  mylwd <- 1
  mycol <- mygrad(length(qs))[i]
  if ((qs[i] == 1) || (qs[i] == 2)) {
    mylwd <- 2.5
    mycol <- "blue"
  }
  lines(loss_lagr_opt_list[[i]] ~ rev(lams), pch = 19, cex = 0.5, col = mycol,
        lwd = mylwd)
}

# auglag
plot(loss_constr_opt_list[[1]] ~ ks, pch = 19, cex = 0.5, col = mygrad(length(qs))[1], type = 'l')
for (i in 2:length(qs)) {
  mylwd <- 1
  mycol <- mygrad(length(qs))[i]
  if ((qs[i] == 1) || (qs[i] == 2)) {
    mylwd <- 2.5
    mycol <- "blue"
  }
  lines(loss_constr_opt_list[[i]] ~ ks, pch = 19, cex = 0.5, col = mycol, lwd = mylwd)
}










