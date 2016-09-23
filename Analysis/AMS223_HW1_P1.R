# AMS 223 Time Series 
# Cheng-Han Yu, Dept of Applied Math and Statistics, UC Santa Cruz
# HW1 Time series and Baysian inference overview
# Problem 1: Chap 1 problem 2

## ----Q1data ----
phi <- 0.9
v <- 1
y0 <- 0.1
n <- 100
y <- rep(NA, n)
# y[1] <- y0
set.seed(123456)
for (i in 1:n) {
    if (i == 1) {
        y[i] <- phi * y0 + rnorm(1, 0, v)
    } else {
        y[i] <- phi * y[i - 1] + rnorm(1, 0, v)
    }
}

## ---- Q1a ----
# (a) MLE for conditional likelihood can be derived analytically
ym1 <- y[-n]
yt <- y[-1]
phi_cmle <- sum(yt * ym1) / (sum(ym1 ^ 2))
v_cmle <- sum((yt - phi_cmle * ym1) ^ 2) / (n - 1)



## ---- Q1b ----
# (b) MLE for the unconditional likelihood using Newton-Raphson method 
# Simulation data from AR(1) with phi = 0.9, v = 1 and y0 = 0.1 
# with sample size n

# Newton-Raphson iteration starting value
phi0 <- 0.8
v0 <- 0.8
theta0 <- c(phi0, v0)
# Hessian(theta0[1], theta0[2])
# solve(Hessian(theta0[1], theta0[2]))
# gradient of the objective function (1.17)
Qstar <- y[1] ^ 2 * (1 - phi ^ 2) + sum((yt - phi * ym1) ^ 2)
gradient <- function(phi, v) {
    dphi <- (-2 * phi / (1 - phi ^ 2)) + (2 / v) * 
        (y[1] ^ 2 * phi + sum(yt * ym1) - phi * sum(ym1 ^ 2))
    dv <- (-n / v) + (1 / v ^ 2) * Qstar
    return (c(dphi, dv))
}

# Hessian of the objective function (1.17)
Hessian <- function(phi, v) {
    dphiphi <- (-2 * (1 + phi ^ 2) / ((1 - phi ^ 2) ^ 2)) + 
        (2 / v) * (y[1] ^ 2 - sum(ym1 ^ 2))
    dvv <- (n / (v ^ 2)) - (2 * Qstar / (v ^ 3))
    dphiv <- (-2 / v ^ 2) * (y[1] ^ 2* phi + sum(yt * ym1) - 
                                phi * sum(ym1 ^ 2))
    return (matrix(c(dphiphi, dphiv, dphiv, dvv), nrow = 2))
}

# Newton-Raphson iteration
theta1 <- theta0 - solve(Hessian(theta0[1], theta0[2])) %*% 
                                gradient(theta0[1], theta0[2])
count <- 1
while (sum(theta1 - theta0) ^ 2 > 1e-8){
    theta0 <- theta1
    theta1 <- theta0 - solve(Hessian(theta0[1], theta0[2])) %*% 
        gradient(theta0[1], theta0[2])
    count <- count + 1
    cat("Iteration = ", count, "\n", sep = "")
    cat("The MLE for (phi, v) = (", theta1[1], ", ", theta1[2], ")", "\n", 
        sep = "")
}
