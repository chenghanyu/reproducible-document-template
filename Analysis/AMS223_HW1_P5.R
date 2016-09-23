# AMS 223 Time Series HW1 Q5 Chap 1 problem 7
# Cheng-Han Yu, Dept of Statistics UC Santa Cruz
# Time series and Baysian inference overview
rm(list = ls())
## ---- Q5a ----
# Metropolis-Hasting algorithm for Threshold autoregressive (TAR) model
# model(1.1) true parameter values
true_phi1 <- 0.9
true_phi2 <- -0.3
true_v <- 1
true_theta <- -1.5

# sample size
n <- 1500

# generate data from model (1.1)
y0 <- 1  # arbitrary initial value
delta <- rep(0, n)  # indicator variable (1 from M1 and 2 from M2)
y <- rep(0, n)
x <- y0
for (i in 1:n) {
    if (x > -true_theta) {
        y[i] <- true_phi1 * x + rnorm(1, 0, true_v)
        delta[i] <- 1
        x <- y[i]
    } else {
        y[i] <- true_phi2 * x + rnorm(1, 0, true_v)
        delta[i] <- 2
        x <- y[i]
    }
}

## ---- Q5a_plot ----
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
plot(y, type = 'l', axes = F, xlab = 'time')
axis(1); axis(2)
plot(delta, type = 'l', axes = F, xlab = 'time')
axis(1); axis(2)

## ---- Q5MCMCalgo ----
# set yt and y_{t-1} vector
yt <- y[-1]
ym1 <- y[-n]

# choose prior parameters
c <- 1
a <- 3
alpha0 <- 3
beta0 <- 0.003

# stotrage
m <- 21000
Theta <- rep(NA, m)
Phi1 <- rep(NA, m)
Phi2 <- rep(NA, m)
V <- rep(NA, m)

# initial values
Phi1[1] <- 0.5
Phi2[1] <- -0.5
Theta[1] <- -2
V[1] <- 2

# counting variable
accept <- 0
count <- 0

Q_fcn <- function(yt, ym1, phi1, phi2, theta) {
    cond <- as.numeric((theta + ym1) > 0)
    sum((yt - ym1 * ifelse(cond, phi1, phi2)) ^ 2)
}

update_phi1 <- function(yt, ym1, phi2, theta, v) {
    cond <- as.numeric((theta + ym1) > 0)
    idx_one <- which(cond == 1)
    sig2_phi1 <- c * v / (c * sum(ym1[idx_one] ^ 2) + v) 
    mu_phi1 <- sig2_phi1 * sum(yt[idx_one] * ym1[idx_one]) / v
    rnorm(1, mu_phi1, sig2_phi1)
}

update_phi2 <- function(yt, ym1, phi1, theta, v) {
    cond <- as.numeric((theta + ym1) > 0)
    idx_zero <- which(cond == 0)
    sig2_phi2 <- c * v / (c * sum(ym1[idx_zero] ^ 2) + v) 
    mu_phi2 <- sig2_phi2 * sum(yt[idx_zero] * ym1[idx_zero]) / v
    rnorm(1, mu_phi2, sig2_phi2)
}

update_v <- function(yt, ym1, phi1, phi2, theta) {
    alpha_v <- (n - 1) / 2 + alpha0
    beta_v <- Q_fcn(yt, ym1, phi1, phi2, theta) / 2 + beta0
    pscl::rigamma(1, alpha_v, beta_v)
}

lpost_theta <- function(yt, ym1, phi1, phi2, theta, v){
    -0.5 * Q_fcn(yt, ym1, phi1, phi2, theta) / v
}

library(pscl)

# phi2 <- Phi2[1]
# theta <- Theta[1]
# v <- V[1]
# for (i in 2:m) { 
#     cat("iter:", i, "\r")
#     
#     # sample phi1 
#     phi1 <- update_phi1(yt, ym1, phi2, theta, v)
#     
#     # sample phi2
#     phi2 <- update_phi2(yt, ym1, phi1, theta, v)
#     
#     # sample v
#     v <- update_v(yt, ym1, phi1, phi2, theta)
#     
#     # use random walk proposal newtheta = theta + N(0, 1) to update theta
#     new.theta <- Theta[i - 1] + rnorm(1, 0, .02)
#     if (new.theta < -a || new.theta > a) {
#         theta <- Theta[i - 1]
#     } else {
#         count <- count + 1
#         u <- runif(1)
#         if (log(u) < (lpost_theta(yt, ym1, phi1, phi2, new.theta, v)
#                      - lpost_theta(yt, ym1, phi1, phi2, Theta[i - 1], v))) {
#             theta <- new.theta
#             accept <- accept + 1
#         }
#     }
#     
#     # store results
#     Phi1[i] <- phi1
#     Phi2[i] <- phi2
#     Theta[i] <- theta
#     V[i] <- v
# }
# draws <- MCMCalgo(yt, ym1, m = 11000)
# burnning and thining
# burn <- 1000
# thin <- 10
# sampleidx = seq(from = (burn + thin), to = m, by = thin)

# trace plots, ACF and histograms
# par(mfrow = c(4, 3), mar = c(4, 4, 4, 1))
# plot(draws[, "Theta"][sampleidx], type = 'l', ylab = "", 
#      main = expression(paste("Trace of ", theta)))
# # abline(h = true_theta, col = 2, lwd = 2)
# acf(draws[, "Theta"][sampleidx], main = expression(paste("ACF of ", theta)))
# hist(draws[,"Theta"][sampleidx], freq = F, breaks = 30, main = "",
#      xlab = expression(theta), col = "navy", border = FALSE)

# main M-H alogorithm
MCMCalgo <- function(yt, ym1, init_phi2, init_theta, init_v, m) {
    # stotrage
    Theta <- rep(NA, m)
    Phi1 <- rep(NA, m)
    Phi2 <- rep(NA, m)
    V <- rep(NA, m)
    
    # counting variable
    accept <- 0
    count <- 0
    
    phi2 <- init_phi2
    theta <- init_theta
    v <- init_v
    
    for (i in 1:m) { 
        # cat("iter:", i, "\r")
        
        # sample phi1 
        phi1 <- update_phi1(yt, ym1, phi2, theta, v)
        
        # sample phi2
        phi2 <- update_phi2(yt, ym1, phi1, theta, v)
        
        # sample v
        v <- update_v(yt, ym1, phi1, phi2, theta)
        
        # use random walk proposal newtheta = theta + N(0, 1) to update theta
        new.theta <- theta + rnorm(1, 0, .025)
        if (new.theta < -a || new.theta > a) {
            theta <- theta
        } else {
            count <- count + 1
            u <- runif(1)
            if (log(u) < (lpost_theta(yt, ym1, phi1, phi2, new.theta, v)
                          - lpost_theta(yt, ym1, phi1, phi2, theta, v))) {
                theta <- new.theta
                accept <- accept + 1
            }
        }
        
        # store results
        Phi1[i] <- phi1
        Phi2[i] <- phi2
        Theta[i] <- theta
        V[i] <- v
    }
    return(list(Phi1 = Phi1, Phi2 = Phi2, V = V, Theta = Theta, 
                accept = accept, count = count))
}

mcmc_sample <- MCMCalgo(yt, ym1, init_phi2 = rnorm(1, 0, 0.25), 
                  init_theta = runif(1, -3, 3), init_v = rigamma(1, 3, 1), 
                  m = 21000) 
# burnning and thining
burn <- 1000
thin <- 10
sampleidx = seq(from = (burn + thin), to = m, by = thin)

post.sample <- function(data, sampleidx) {
    Phi1 <- data$Phi1[sampleidx]
    Phi2 <- data$Phi2[sampleidx]
    V <- data$V[sampleidx]
    Theta <- data$Theta[sampleidx]
    draws <- cbind(Phi1, Phi2, V, Theta)
    colnames(draws) <- c("phi_1", "phi_2", "v", "theta")
    return(draws)
}

draws <- post.sample(mcmc_sample, sampleidx)

## ---- Q5MCMCplot ----
par(mfrow = c(4, 3), mar = c(4, 4, 4, 1))
plot(draws[, "phi_1"], type = 'l', ylab = "", 
     main = expression(paste("Trace of ", phi[1])))
# abline(h = true_phi1, col = 2, lwd = 2)
acf(draws[, "phi_1"], main = expression(paste("ACF of ", phi[1])))
hist(draws[, "phi_1"], freq = F, breaks = 30, main = "",
     xlab = expression(phi[1]), col = "navy", border = FALSE)

plot(draws[, "phi_2"], type = 'l', ylab = "", 
     main = expression(paste("Trace of ", phi[2])))
# abline(h = true_phi2, col = 2, lwd = 2)
acf(draws[, "phi_2"], main = expression(paste("ACF of ", phi[2])))
hist(draws[, "phi_2"], freq = F, breaks = 30, main = "",
     xlab = expression(phi[2]), col = "navy", border = FALSE)

plot(draws[, "v"], type = 'l', ylab = "", 
     main = expression(paste("Trace of ", v)))
# abline(h = true_v, col = 2, lwd = 2)
acf(draws[, "v"], main = expression(paste("ACF of ", v)))
hist(draws[, "v"], freq = F, breaks = 30, main = "",
     xlab = expression(v), col = "navy", border = FALSE)

plot(draws[, "theta"], type = 'l', ylab = "", 
     main = expression(paste("Trace of ", theta)))
# abline(h = true_theta, col = 2, lwd = 2)
acf(draws[, "theta"], main = expression(paste("ACF of ", theta)))
hist(draws[, "theta"], freq = F, breaks = 30, main = "",
     xlab = expression(theta), col = "navy", border = FALSE)

## ---- Q5MCMCresult ----
# acceptance rate
accept_rate = mcmc_sample$accept / mcmc_sample$count

quan025 <- function(x) {
    quantile(x, prob = 0.025)
}
quan975 <- function(x) {
    quantile(x, prob = 0.975)
}

library(coda)
# draws <- cbind(Phi1, Phi2, V, Theta)[sampleidx, ]
colnames(draws) <- c("$\\phi_1$", "$\\phi_2$", "$v$", "$\\theta$")
result <- round(cbind(apply(draws, 2, effectiveSize), 
                      apply(draws, 2, mean), 
                      apply(draws, 2, quan025), 
                      apply(draws, 2, quan975)), 4)
colnames(result) <- c("effective size", "post. mean", "2.5% quantile", 
                      "97.5% quantile")
library(xtable)
print(xtable(result, caption = "Posterior Summary", label = "MCMCresult"),
      sanitize.rownames.function=function(x){x})
# for html output
# print(xtable(result, caption = "Posterior Summary", label = "MCMCresult"),
#       sanitize.rownames.function=function(x){x}, type = "html")

## ---- Q5coda ----
library(coda)
library(parallel)
library(doParallel)
detectCores()
cl <- makeCluster(2, type = "FORK")
registerDoParallel(cl)
getDoParWorkers()

MCMCalgo.mc <- function(s) {
    # initial values
    init_phi2 <- rnorm(5, 0, 0.25)
    init_theta <- runif(5, -3, 3)
    init_v <- rigamma(5, 3, 1)
    MCMCalgo(yt, ym1, init_phi2[s], init_theta[s], init_v[s], m = 21000)
}

system.time(draws_mc <- mclapply(1:5, MCMCalgo.mc, mc.cores = 2))

stopCluster(cl)

#------------------------------------------------------------------------------
# Analysis using coda package
#------------------------------------------------------------------------------


draws.mc = lapply(draws_mc, post.sample, sampleidx = sampleidx)
coda.draws.mc = lapply(draws.mc, mcmc)
# mean, sd, and quantiles of the two chains
lapply(coda.draws.mc, summary)
# traceplots, and densities 
lapply(coda.draws.mc, plot)
# pairwise correlations
lapply(coda.draws.mc, function(x) pairs(data.frame(x)))

# Convergence diagnostics
#------------------------
combinedchains = mcmc.list(coda.draws.mc[[1]], coda.draws.mc[[2]], 
                           coda.draws.mc[[3]], coda.draws.mc[[4]],
                           coda.draws.mc[[5]])
plot(combinedchains)
# acf
autocorr.diag(combinedchains)
autocorr.plot(combinedchains)
# crosscorr
crosscorr.plot(combinedchains)
# Gelman and Rubin potential scale reduction factor 
gelman.diag(combinedchains)  # should be close to 1
gelman.plot(combinedchains)
# Geweke’s convergence diagnostic
geweke.diag(combinedchains)  # Z-scores
geweke.plot(combinedchains)
# Heidelberger and Welch’s convergence diagnostic
heidel.diag(combinedchains)

# Raftery and Lewis’s diagnostic
raftery.diag(combinedchains)


