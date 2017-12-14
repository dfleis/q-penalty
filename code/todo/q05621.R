##### setup #####
library(alabama)
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
beta1 <- 0.75
beta2 <- 0.95
b1min <- 0
b1max <- 1
nb1 <- 1e2
b2min <- 0
b2max <- 1
nb2 <- 1e2

kmin <- 0
kmax <- 3
nks <- 200

q <- 0.56214641035#01893229543725289214607

lammin <- 0
lammax <- 2
nlams <- 1e4

lams1 <- seq(lammin, 1/2, length.out = round(0.05 * nlams))
lams2 <- seq(1/2, lammax, length.out = round(0.95 * nlams))
lams1 <- lams1[which(lams1 != 1/2)]
lams2 <- lams2[which(lams2 != 1/2)]
lams <- c(lams1, lams2)

b1 <- seq(b1min, b1max, length.out = nb1)
b2 <- seq(b2min, b2max, length.out = nb2)
ks <- seq(kmin, kmax, length.out = nks)
qs <- c(q, 1, 2)

xmin <- b1min
xmax <- b1max * 1.3
ymin <- b2min
ymax <- b2max * 1.3
b1plot <- seq(xmin, xmax, length.out = nb1)
b2plot <- seq(ymin, ymax, length.out = nb2)

##### define some functions #####
obj_fxn <- function(x, y) -(x * beta1 + y * beta2) + 1/2 * (x^2 + y^2)
obj_grad_fxn <- function(x, y) c(-beta1 + x, -beta2 + y)
constr_fxn <- function(x, y, k, q) k - abs(x)^q - abs(y)^q # hin(x) >= 0
constr_fxn2 <- function(x, y, k, q) abs(x)^q + abs(y)^q - k # heq(x) == 0
constr_fxn3 <- function(x, y, k) (x - k)^2 + (y - k)^2 - k^2 # heq(x) == 0


sol2x <- function(a, b, l) {
  (-2 * b * l + 2 * sqrt(2) * (-1 + 2 * l) * 
     sqrt((l^2 * (a * b + 2 * a^2 * (-1 + l) * l + 2 * b^2 * (-1 + l) * l))/(1 - 2 * l)^2) +
     a * (1 + 2 * l - 4 * l^2))/(1 + 2 * l - 12 * l^2 + 8 * l^3)
}
sol2y <- function(a, b, l) {
  (-2 * a * l + 2 * sqrt(2) * (-1 + 2 * l) * 
     sqrt((l^2 * (b * a + 2 * b^2 * (-1 + l) * l + 2 * a^2 * (-1 + l) * l))/(1 - 2 * l)^2) +
     b * (1 + 2 * l - 4 * l^2))/(1 + 2 * l - 12 * l^2 + 8 * l^3)
}

##### compute #####
loss_baseline <- outer(b1plot, b2plot, function(x, y) obj_fxn(x, y))
loss_opt_list <- list()

# auglag: augmented lagragian algorithm
pt <- proc.time()
for (q_idx in 1:length(qs)) {
  loss_opt_ks <- sapply(ks, function(k)
    alabama::auglag(par = c(0, 0),
                    fn  = function(params) obj_fxn(params[1], params[2]),
                    gr  = function(params) obj_grad_fxn(params[1], params[2]),
                    hin = function(params) constr_fxn(params[1], params[2], k, qs[q_idx]),
                    control.outer = list(maxeps = 1e-4, itmax = 1e4, trace = F))$par)
  loss_opt_list[[q_idx]] <- loss_opt_ks
}
proc.time() - pt

pt <- proc.time()
loss_opt_circ <- sapply(ks, function(k) 
  alabama::auglag(par = c(0, 0),
                  fn  = function(params) obj_fxn(params[1], params[2]),
                  gr  = function(params) obj_grad_fxn(params[1], params[2]),
                  hin = function(params) constr_fxn3(params[1], params[2], k),
                  control.outer = list(maxeps = 1e-4, itmax = 1e4, trace = F))$par)
proc.time() - pt

# remove infeasible points
loss_opt_circ[2,][which((loss_opt_circ[2,] > beta2) & (loss_opt_circ[1,] < beta1))] <- 0
loss_opt_circ[2,][which((loss_opt_circ[2,] > beta2) & (loss_opt_circ[1,] > beta1))] <- beta2

loss_opt_circ[1,][which(loss_opt_circ[1,] > beta1)] <- beta1
loss_opt_circ[1,][which((loss_opt_circ[1,] > 0) & (loss_opt_circ[2,] == 0))] <- 0

# compute solution
s2x <- sol2x(beta1, beta2, lams)
s2y <- sol2y(beta1, beta2, lams)

# remove infeasible points

# s2y[which((s2y > beta2) & (s2x < beta1))] <- 0
# s2y[which((s2y > beta2) & (s2x < beta1))] <- beta2
# 
# s2x[which((s2x > beta1)] <- beta1
# s2x[which((s2x > 0) & (s2y == 0))] <- 0
# 
# s2x[which(s2x - s2y > 0)] <- NA # remove bottom branch

##### plots #####
contour(b1plot, b2plot, loss_baseline, 
        lwd = 0.75, nlevels = 20, drawlabels = F, col = 'gray50',
        xlim = c(xmin, xmax), 
        ylim = c(ymin, ymax),
        xlab = expression(beta[1]),
        ylab = expression(beta[2]),
        main = substitute(
          paste("q = ", q),
          list(q = q)))
lines(
  x = loss_opt_circ[1, order(loss_opt_circ[1,])],
  y = loss_opt_circ[2, order(loss_opt_circ[1,])],
  lwd = 2.5, col = 'blue')
for (i in 1:length(qs)) {
  mylwd <- 1
  mycol <- "red"
  mylty <- 'solid'
  if ((qs[i] == 1) || (qs[i] == 2)) {
    mycol <- "black"
    mylwd <- 1.5
    if (qs[i] == 1) 
      mylty <- 'dotdash'
    if (qs[i] == 2)
      mylty <- 'longdash'
  }
  lines(loss_opt_list[[i]][2,] ~ loss_opt_list[[i]][1,], col = mycol,
        type = 'l', pch = 19, lwd = mylwd, lty = mylty)
}
points(beta1, beta2, pch = 19, cex = 0.75, col = 'red')
text(beta1 + 0.09, beta2 + 0.02, labels = expression(beta^OLS))
lines(s2x, s2y, 
      col = 'green3', lty = 'dashed', lwd = 2)
legend("bottomright", 
       legend = c("q norm", "q norm (cncv aprx)", "q norm (cncv aprx) soln", "lasso", "ridge"),
       col = c("red", "blue", "green3", "black", "black"), lwd = c(1, 2, 2, 1.5, 1.5), 
       seg.len = 2, lty = c("solid", "solid", "dashed", "dotdash", "longdash"))




