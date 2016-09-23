################################################################################
# AMS 223 Time Series 
# Cheng-Han Yu, Dept of Applied Math and Statistics, UC Santa Cruz
# HW1 Time series and Baysian inference overview
# Problem 4: Chap 1 problem 5
################################################################################
## ---- Q4ts ----
# (a) Sample 200 obs
n <- 200
v <- 1
# choose model 1 parameters
phi1 <- 0.5
phi2 <- 0.15
# (a.1) sample 200 obs from model 1
# set.seed(1234)
y <- stats::arima.sim(n = n, model = list(ar = c(phi1, phi2)), sd = v)
# head(y)
plot(y, type = "l", xlab = "time", ylab = "", axes = F)
axis(1)
axis(2)

a <- 1
b <- 2
w0 <- 0.3
# (a.2) sample 200 obs from model 2
angle = 2 * w0 * pi * 1:n
x = a * cos(angle) + b * sin(angle) + rnorm(n, 0, v)
# head(x)
plot(x, type = "l", xlab = "time", ylab = "", axes = F)
axis(1)
axis(2)


## ---- Q4b1 ----
# (b) Find MLE
# (b.1) MLE for model 1
f1 <- y[-c(1, n)]  # 1st col of F'
f2 <- y[-c(n - 1, n)]  # 2nd col of F'
Ft <- cbind(f1, f2)  # Ft = F' in the model
mle1 <- chol2inv(chol(crossprod(Ft))) %*% t(Ft) %*% y[-c(1, 2)]
R1 <- sum((y[3:n] - Ft %*% mle1) ^ 2)
mle_v1 <- R1 / (n - 2)  # sample size n-2
s2_1 <- R1 / (n - 2 - 2)  # 2 parameters phi1, phi2


## ---- Q4b2 ----
# (b.2) MLE for model 2 
# depend on Q4a2
Xt <- cbind(cos(angle), sin(angle))
mle2 <- chol2inv(chol(crossprod(Xt))) %*% t(Xt) %*% x

R2 <- sum((x - Xt %*% mle2) ^ 2)
mle_v2 <- R2 / (n)  # sample size n

s2_2 <- R2 / (n - 2) # 2 parameters a, b



## ---- Q4c1 ----
# (c) Find MAP
# (c.1) MAP for model 1
# (phi1, phi2) same as MLE
map1 <- mle1

# v is IG((n - 2 - 2) / 2), (n - 2 - 2) * s ^ 2 / 2).
# v_map = mode of IG = ((n - 2 - 2) * s ^ 2 / 2) / ((n - 2 - 2) / 2 + 1)
map_v1 <- ((n - 2 - 2) * s2_1 / 2) / ((n - 2 - 2) / 2 + 1)


## ---- Q4c2 ----
# (c.2) MAP for model 2
# (phi1, phi2) same as MLE
map2 <- mle2

# v is IG((n - 2) / 2), (n - 2) * s ^ 2 / 2).
# v_map = (n - 2) * s ^ 2 / 2 /((n - 2) / 2 + 1)
map_v2 <- ((n - 2) * s2_2 / 2) / ((n - 2) / 2 + 1)

## ---- Q4d ----
suppressMessages(library(pscl))
library(tikzDevice)
alpha1 <- (n - 4) / 2
beta1 <- (n - 4) * s2_1/2
v1 <- seq(0.1, 3, 0.01)
sample1 <- pscl::rigamma(500, alpha1, beta1)
hist_v <- function(samp, mdl, prior.type, alpha, beta, mle, s2, map) {
    par(mar = c(4, 4, 2, .1))
    hist(samp, prob = T, 
         main = paste("Marginal posterior of $v$: model", mdl, prior.type),
         xlab = '$v$', ylab = '$p(v|\\mathbf{y})$', breaks = 50, 
         col  ="lightgray", cex.lab = 1.5, cex.main = 1.5, 
         border = "white")
    lines(v1, densigamma(v1, alpha, beta), type = 'l')
    points(mle, 0, pch = 16, col = "red")  # MLE
    points(s2, 0.1, pch = 17, col = "blue")   # s^2
    points(map, 0.2, pch = 15, col = "brown")  # MAP
    legend("topright", title = "Estimator", bty = "n",
           c("MLE", expression(s^2), 
             substitute(paste("MAP"[prior.type]))), 
           pch = c(16, 17, 15), col = c("red", "blue", "brown"))
}
hist_v(sample1, 1, "ref", alpha1, beta1, mle_v1, s2_1, map_v1)

library(mvtnorm)
library(fields)
den_coef <- function(k, map, Sigma, mdl, prior.type) {
    m <- length(k)
    mu <- as.vector(map) 
    # Omega <- s2 * chol2inv(chol(crossprod(design_mat)))
    Z = matrix(NA, nrow = m, ncol = m)
    for (i in 1:m) {
        for (j in 1:m) {
            Z[i, j] = dmvt(c(k[i], k[j]), delta = mu, sigma = Sigma, df = n - 4, 
                           log = FALSE)
        }
    }
    # Z = outer(k, k, dmvt, delta = mu, sigma = Sigma, df = n - 4, log = FALSE) 
    # does not work. Why?
    par(mar = rep(.05, 4))
    persp(k, k, Z, theta = -40, phi = 30, col = "lightblue", border = "darkgray",
          box = T, axes = F) 
    par(mar = rep(3, 4))
    image.plot(k, k, Z, axes = T, xlab = "", ylab = "", cex.axis = 0.8)
    contour(k, k, Z, add = TRUE, col = "white")
    title(main = paste("Student-t for", 
                       ifelse(mdl == 1, "$(\\phi_1, \\phi_2)$", "(a, b)"), 
                       "\n model", mdl, prior.type), cex.main = 1)
}
k <- seq(-0.2, 0.8, length = 50)
Omega <- s2_1 * chol2inv(chol(crossprod(Ft)))
par(mfrow = c(1, 2))
den_coef(k, mle1, Omega, 1, "ref")
# dev.off()

## ---- Q4e2v ----
# (e) model 2: sketch marginal posterior of v and (a, b)
# sketch v
alpha2 <- (n - 2) / 2
beta2 <- (n - 2) * s2_2 / 2
sample2 <- rigamma(500, alpha2, beta2)
hist_v(sample2, 2, "ref", alpha2, beta2, mle_v2, s2_2, map_v2)


## ---- Q4e2phi ----
# sketch (a, b)
k2 <- seq(0.6, 2.4, length = 50)
Omega2 <- s2_2 * chol2inv(chol(crossprod(Xt)))
den_coef(k2, mle2, Omega2, 2, "ref")


## ---- Q4e2 ----
# (e) model 2: sketch marginal posterior of v and (a, b)
# sketch v
alpha2 <- (n - 2) / 2
beta2 <- (n - 2) * s2_2 / 2
sample2 <- rigamma(500, alpha2, beta2)
hist_v(sample2, 2, "ref", alpha2, beta2, mle_v2, s2_2, map_v2)

# sketch (a, b)
k2 <- seq(0.6, 2.4, length = 50)
Omega2 <- s2_2 * chol2inv(chol(crossprod(Xt)))
par(mfrow = c(1, 2))
den_coef(k2, mle2, Omega2, 2, "ref")


## ---- Q4f1c ----
# (f) conjugate prior case: redo (c), (d), and (e)
# model 1
# set up conjugate priors
# (phi1, phi2) ~ N(m0, vC0), v ~ IG(n0/2, d0/2)
m01 <- c(0.2, -0.5)
C01 <- diag(3, 2)
n01 <- 10
d01 <- 20

# redo(c) posterior of (phi1, phi2) ~ T(m, vC) df = nstar 
# and v ~ IG(nstar/2, dstar/2)
Q1 <- Ft %*% C01 %*% t(Ft) + diag(1, n - 2)
e1 <- y[3:n] - Ft %*% m01
AA <- C01 %*% t(Ft) %*% chol2inv(chol(Q1))
m1 <- m01 + AA %*% e1
C1 <- C01 - AA %*% Ft %*% C01

nstar1 <- (n - 2) + n01
dstar1 <- t(e1) %*% chol2inv(chol(Q1)) %*% e1 + d01

map1_conj <- m1
map_v1_conj <- (dstar1 / 2) / ((nstar1 / 2) + 1)


## ---- Q4f1d ----
# redo(d) sketch marginal posterior of v and (phi1, phi2)
# sketch v
sample1conj <- rigamma(500, nstar1 / 2, dstar1 / 2)
hist_v(sample1conj, 1, "conj", nstar1 / 2, dstar1 / 2, mle_v1, 
       s2_1, map_v1_conj)

# sketch (phi1, phi2)
Omega <- as.vector(dstar1 / nstar1) * C1 
par(mfrow = c(1, 2))
den_coef(k, map1_conj, Omega, 1, "conj")

## ---- Q4f2c ----
# model 2
# set up conjugate priors
# (a, b) ~ N(m0, vC0), v ~ IG(n0/2, d0/2)
m02 <- c(2, 3)
C02 <- diag(1, 2)
n02 <- 10
d02 <- 20

# redo(c) posterior of (phi1, phi2) ~ T(m, vC) df = nstar 
# and v ~ IG(nstar/2, dstar/2)
Q2 <- Xt %*% C02 %*% t(Xt) + diag(1, n)
e2 <- (x - Xt %*% m02)
BB <- C02 %*% t(Xt) %*% chol2inv(chol(Q2))
m2 <- m02 + BB %*% e2
C2 <- C02 - BB %*% Xt %*% C02

nstar2 <- n + n02
dstar2 <- t(e2) %*% chol2inv(chol(Q2)) %*% e2 + d02

map2_conj <- m2
map_v2_conj <- (dstar2 / 2) / ((nstar2 / 2) + 1)


## ---- Q4f2e ----
# redo(e) sketch marginal posterior of v and (a, b)
# sketch v
v2 <- seq(0.1, 3, 0.01)
sample2conj <- rigamma(500, nstar2/2, dstar2/2)
hist_v(sample2conj, 2, "conj", nstar2 / 2, dstar2 / 2, 
       mle_v2, s2_2, map_v2_conj)

# sketch (a, b)
Omega2 <- as.vector(dstar2 / nstar2) * C2
par(mfrow = c(1, 2))
den_coef(k2, map2_conj, Omega2, 2, "conj")


## ---- Q4f1senv ----
# Sensitivity analysis
# model 1 
# reset conjugate priors
m01new <- c(0.9, -0.1)
C01new <- diag(2, 2)
n01new <- 50 
d01new <- 100

# redo(c) 
Q1new <- Ft %*% C01new %*% t(Ft) + diag(1, n-2)
e1new <- (y[3:n] - Ft %*% m01new)
AAnew <- C01new %*% t(Ft) %*% chol2inv(chol(Q1new))
m1new <- m01new + AAnew %*% e1new
C1new <- C01new - AAnew %*% Ft %*% C01new

nstar1new <- (n - 2) + n01new
dstar1new <- t(e1new) %*% chol2inv(chol(Q1new)) %*% e1new + d01new

map1_conj_new <- m1new
# compare map
# map1_conj
# map1_conj_new
map_v1_conj_new <- (dstar1new / 2) / ((nstar1new / 2) + 1)
# map_v1_conj
# map_v1_conj_new
par(mfrow = c(1, 1), mar = c(4, 4, 1, 1))
plot(v1, densigamma(v1, nstar1 / 2, dstar1 / 2), type = 'l', 
     xlab = '$v$', ylab = '$p(v|\\mathbf{y})$', lwd = 2, col = 2)
lines(v1, densigamma(v1, nstar1new / 2, dstar1new / 2), lwd = 3, col = 4)
legend("topright", title = "Sensitivity of model 1: $v$", 
       c("$n_0 = 10$, $d_0 = 20$", "$n_0 = 50$, $d_0 = 100$"), lwd = c(2, 3), 
       cex = 0.9, bty = "n", col = c(2, 4))

## ---- Q4f1senphi ----
Omeganew <- as.vector(dstar1new / nstar1new) * C1new
par(mfrow = c(2, 2))
den_coef(k, map1_conj, Omega, 1, "$\\mathbf{m_0} = (0.2, -0.5)$")
den_coef(k, map1_conj_new, Omega, 1, "$\\mathbf{m_0} = (0.9, -0.1)$")

## ---- Q4f2senv ----
# model 2
# reset conjugate priors
m02new <- c(2, 3) # = m02
m02new1 <- c(10, 20)
C02new <- diag(1, 2)
n02new <- 50
d02new <- 100

# redo(c) 
Q2new <- Xt %*% C02new %*% t(Xt) + diag(1, n)
e2new <- (x - Xt %*% m02new)
e2new1 <- (x - Xt %*% m02new1)
BBnew <- C02new %*% t(Xt) %*% chol2inv(chol(Q2new))
m2new <- m02new + BBnew %*% e2new
m2new1 <- m02new1 + BBnew %*% e2new1
C2new <- C02new - C02new %*% t(Xt) %*% chol2inv(chol(Q2new)) %*% Xt %*% C02new

nstar2new <- (n) + n02new
nstar2new1 <- (n) + n02new

dstar2new <- t(e2new) %*% chol2inv(chol(Q2new)) %*% e2new + d02new
dstar2new1 <- t(e2new1) %*% chol2inv(chol(Q2new)) %*% e2new1 + d02new

map2_conj_new1 <- m2new1
map_v2_conj_new <- (dstar2new / 2) / ((nstar2new / 2) + 1)


par(mfrow = c(1, 1), mar = c(4, 4, 1, 1))
plot(v2, densigamma(v2, nstar2 / 2, dstar2 / 2), type = 'l', 
     xlab = '$v$', ylab = '$p(v|\\mathbf{y})$', lwd = 2, col = 2)
lines(v2, densigamma(v2, nstar2new/2, dstar2new/2), lwd = 3, col = 4)
legend("topright", bty = "n",
       title = "Sensitivity of model 2: $v$ and same $\\mathbf{m_0} = (2, 3)'$",
       c("$n_0 = 10$, $d_0 = 20$", "$n_0 = 50$, $d_0 = 100$"), lwd = c(2, 3), 
       cex = 0.9, col = c(2, 4))

## ---- Q4f2senab1----
mu <- as.vector(m2new)
munew <- as.vector(m2new1)
Omeganew <- as.vector(dstar2new / nstar2new) * C2
Omeganew1 <- as.vector(dstar2new1 / nstar2new1) * C2new

library(fMultivar)
Z = matrix(dmvst(K2, 2, mu, Omega, df = nstar2), length(k2))
Znew = matrix(dmvst(K2, 2, munew, Omeganew, df = nstar2new), length(k2))
par(mfrow = c(2, 2))
den_coef(k2, m2new, Omeganew1, 2, "same $(n_0, d_0) = (50, 100)$")
den_coef(k2, m2new1, Omeganew1, 2, "same $(n_0, d_0) = (50, 100)$")

## ---- Q4f2senab2----
par(mfrow = c(2, 2))
den_coef(k2, m2new, Omeganew, 2, "$(n_0, d_0) = (10, 20)$")
den_coef(k2, m2new1, Omeganew1, 2, "$(n_0, d_0) = (50, 100)$")
