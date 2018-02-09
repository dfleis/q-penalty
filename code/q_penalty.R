###
#
# Comparison of the constrained optimization problem (given p predictors)
# argmin ( (Y - X beta)^T (Y - X beta) ) subject to sum^p_{j = 1} abs(beta_j)^q <= k^q, k > 0
# verus the Lagrangian
# argmin ( (Y - X beta)^T (Y - X beta) + lambda * sum^p_{j = 1} abs(beta_j)^q ), lambda > 0
#
# For simulation purposes we assume the data is orthogonal X^T X = I_p so that
# beta_ols = X^T Y, yielding the equivalent loss functions:
# argmin ( -beta^T beta_ols + beta^T beta ) subject to sum^p_{j = 1} abs(beta_j)^q <= k^q, k > 0
# verus the Lagrangian
# argmin ( -beta^T beta_ols + beta^T beta + lambda * sum^p_{j = 1} abs(beta_j)^q ), lambda > 0
###
library(alabama)
mygrad <- colorRampPalette(c("red", "black"))
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

##### define some parameters ######
# q values
qs <- 0.1 * (1:20)
qs <- c(qs, 1, 2, 3, 4, 5, 10, 20) 

# theoretical OLS estimates
beta_ols1 <- 0.7
beta_ols2 <- 0.9

# domain to compute regularized estimates according to the q-penalty
b1min <- 0
b1max <- 1
db1   <- 0.01
b2min <- 0
b2max <- 1
db2   <- 0.01

# values of k and lambda to compute our regularized estimates over
kmin <- 0
kmax <- 3.5
dk   <- 0.01
lammin <- 0
lammax <- 20
dlam   <- 0.01

qs <- qs[which(qs != 0)]
qs <- sort(unique(qs)) # remove duplicates
b1   <- seq(b1min, b1max, db1)
b2   <- seq(b2min, b2max, db2)
ks   <- seq(kmin, kmax, dk)
lams <- seq(lammin, lammax, dlam)

##### functions ######
obj_fxn <- function(x, y, beta1, beta2) -(x * beta1 + y * beta2) + 1/2 * (x^2 + y^2)
obj_grad_fxn <- function(x, y, beta1, beta2) c(-beta1 + x, -beta2 + y)
constr_fxn <- function(x, y, k, q) k - abs(x)^q - abs(y)^q # hin(x) >= 0
constr_fxn2 <- function(x, y, k, q) abs(x)^q + abs(y)^q - k # g_ineq(x) <= 0
obj_fxn_lagr <- function(x, y, lam, q, beta1, beta2) obj_fxn(x, y, beta1, beta2) + lam * (abs(x)^q + abs(y)^q)

##### DO WORK ######
loss_constr_opt_list <- list()
loss_lagr_opt_list <- list()

# baseline loss contours/surface
loss_ols <- outer(b1, b2, function(x, y) obj_fxn(x, y, beta_ols1, beta_ols2))

pt <- proc.time()
for (q_idx in 1:length(qs)) {
  q <- qs[q_idx]
  
  ## lagragian multiplier framework
  pt <- proc.time()
  loss_lagr_opt_lams <- sapply(lams, FUN = function(l) 
    optim(c(0, 0), # start search at (0, 0)
          fn = function(params) obj_fxn_lagr(params[1], params[2], l, q, beta_ols1, beta_ols2))$par)
  
  ## constrained framework
  # auglag: augmented lagragian algorithm
  loss_constr_opt_ks <- sapply(ks, function(k)
    alabama::auglag(par = c(0, 0), # start searach at (0, 0)
                    fn  = function(params) obj_fxn(params[1], params[2], beta_ols1, beta_ols2),
                    gr  = function(params) obj_grad_fxn(params[1], params[2], beta_ols1, beta_ols2),
                    hin = function(params) constr_fxn(params[1], params[2], k, q),
                    control.outer = list(maxeps = 1e-4, itmax = 1e4, trace = F))$par)
  tm <- proc.time() - pt
  
  loss_lagr_opt_list[[q_idx]] <- loss_lagr_opt_lams
  loss_constr_opt_list[[q_idx]] <- loss_constr_opt_ks
  
  print(paste0("q index = ", round(q_idx/length(qs) * 100, 3), "% ", round(unname(tm[3]), 3), "s"))
}
proc.time() - pt

##### plots #####
## lagragian multiplier
#image(b1, b2, loss_ols, col = parula(1e2),
#      xlab = expression(beta[1]),
#      ylab = expression(beta[2]),
#      main = paste("lagr opt"))
pdf("./img/lagr_opt.pdf", height = 5, width = 6)
contour(b1, b2, loss_ols, col = 'gray50', drawlabel = F,
        xlab = expression(beta[1]),
        ylab = expression(beta[2]),
        main = substitute(paste("Lagrangian Mult., ", 
                                lambda %in% "[", lmin, ", ", lmax, "], ", 
                                q %in% "[", qmin, ", ", qmax, "]"),
                          list(lmin = min(lams), lmax = max(lams), qmin = min(qs), qmax = max(qs))))
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
points(beta_ols1, beta_ols2, pch = 19)
text(beta_ols1 + 0.05, beta_ols2 + 0.075, labels = expression(beta^OLS))
dev.off()

## auglag
#image(b1, b2, loss_ols, col = parula(1e2),
#      xlab = expression(beta[1]),
#      ylab = expression(beta[2]),
#      main = paste("constr opt: auglag"))
pdf("./img/constr_opt.pdf", height = 5, width = 6)
contour(b1, b2, loss_ols, col = 'gray50', drawlabel = F,
        xlab = expression(beta[1]),
        ylab = expression(beta[2]),
        main = substitute(paste("Constrained Optim., ", 
                                k %in% "[", kmin, ", ", kmax, "], ",
                                q %in% "[", qmin, ", ", qmax, "]"),
                          list(kmin = min(ks), kmax = max(ks), qmin = min(qs), qmax = max(qs))))
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
points(beta_ols1, beta_ols2, pch = 19)
text(beta_ols1 + 0.05, beta_ols2 + 0.075, labels = expression(beta^OLS))
dev.off()













