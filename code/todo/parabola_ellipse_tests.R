


halfnorm <- function(x, r) {
  r^2 - 2 * r * sqrt(x) + x
}
ellipse <- function(x, a, b, h, k) {
  #-(1 - (x - h)^2)^(1/2) + k
  k + (b - b/a * (x - h)^2)^(1/2)
}

a <- 1
b <- 1
h <- 1
k <- -1
r <- 1
xellipse <- seq(h - a, h + a, length.out = 1e3)
x <- seq(0.01, 1, length.out = 1e3)
kvals <- -r^2 + 2 * r * sqrt(x) - x - sqrt(1/a * (-b * (-a + h^2 - 2 * h * x + x^2) ))



plot(halfnorm(x, r) ~ x, type = 'l', xlim = c(0, 2), ylim = c(0, 2), lwd = 2)
lines(-ellipse(xellipse, a, b, h, k) ~ xellipse, col = 'red')
lines(ellipse(xellipse, a, b, h, -k) ~ xellipse, col = 'red')
points(h, -k, pch = 19)









