# sumsquare <- function (y, theta, phi1, phi2) {
#   A <- 0
#   y0 <- 6
#   if (y0 > -theta) {
#     A <- (y[1] - phi1*y0)^2
#   } else {
#     A <- (y[1] - phi2*y0)^2
#   }
#   for (i in 1: (M-1)) {
#     if (y[i] > -theta) {
#       A <- A + (y[i+1] - phi1*y[i])^2
#     } else {
#       A <- A + (y[i+1] - phi2*y[i])^2
#     }
#   }
#   return (A)
# }
# index <- function(y, theta) {
#   delta <- rep(NA, 1000) 
#   for (i in 1: 1000) {
#     if (y0 > -theta) {
#       delta[i] <- 1
#       y0 <- y[i]
#     } else {
#       delta[i] <- 2
#       y0 <- y[i]
#     }
#   }
#   idx1 <- which(delta %in% 1)
#   idx2 <- which(delta %in% 2)
#   n1 <- length(idx1)
#   n2 <- length(idx2)
#   previdx1 <- idx1 - rep(1, n1)
#   previdx2 <- idx2 - rep(1, n2)  
#   ymodel1 <- rep(NA, n1) 
#   for (i in ind1) {  
#     ymodel1[i] <- y[i]
#   }
#   ymodel2 <- rep(NA, n2)
#   for (i in ind2) {
#     ymodel2[i] <- y[i]
#   }
#   A <- 0
#   if (y0 > -theta) {
#     A <- (y[1] - phi1*y0)^2
#   } else {
#     A <- (y[1] - phi2*y0)^2
#   }
#   A1 <- A + sum((ymodel1[2: n1] - phi1*ymodel1[1: n-1])^2)
#   A2 <- A + sum((ymodel2[2: n1] - phi2*ymodel2[1: n-1])^2)
#   return(c(A1, A2))
# }

# prevymodel1 <- c(y0, y[prevind1[2:n1]])


# # posterior distribution
# logpost <- function(y, theta, phi1, phi2, v) {
#   cond1 <- as.numeric((theta + ytm1) > 0)
#   cond2 <- as.numeric(!cond1)
#   a1 <- 1/sqrt(v)*exp(-.5/v*(yt-ytm1*phi1)^2)
#   a2 <- 1/sqrt(v)*exp(-.5/v*(yt-ytm1*phi2)^2)
#   #A <- sumsquare(y, theta, phi1, phi2)
#   pst <- v^(-((/2)+alpha0+1))*exp(-(1/(2*v))*A)*exp((-1/2)*phi1^2)*exp((-1/2)*phi1^2)*exp(-beta0/v)
#   return(pst)
# }

# condi.v <- function(y, theta, phi1, phi2) {
#   A <- sumsquare(y, theta, phi1, phi2)
#   alpha0 <- 2
#   beta0 <- 4
#   sim.v <- rigamma(1, alpha=(M/2)+alpha0, beta=(1/2)*A+beta0)
#   return(sim.v)
# }

# condi.phi1 <- function(y, theta, phi2, v) {
#   y0 <- 6
#   x <- y0
#   delta <- rep(NA, M) 
#   for (i in 1: M) {
#     if (x > -theta) {
#       delta[i] <- 1
#       x <- y[i]
#     } else {
#       delta[i] <- 2
#       x <- y[i]
#     }
#   }
#   idx1 <- which(delta %in% 1)
#   n1 <- length(idx1)
#   previdx1 <- idx1 - rep(1, n1)
#   ymodel1 <- rep(NA, n1) 
#   for (i in 1: n1) {  
#     ymodel1[i] <- y[idx1[i]]
#   }
#   C <- 0
#   if (y0 > -theta) {
#     C <- ymodel1[1]*y0
#     previdx1 <- previdx1[2: n1]
#     C1 <- (C + sum(ymodel1[2:n1]*y[previdx1]))/(sum(y[previdx1]^2)+y0+v)
#     sim.phi1 <- rnorm(1, mean = C1, sd = sqrt((sum(y[previdx1]^2)+y0 + v)))
#   } else {
#     C1 <- sum(ymodel1*y[previdx1])/(sum(y[previdx1]^2)+v)
#     sim.phi1 <- rnorm(1, mean = C1, sd = sqrt((sum(y[previdx1]^2)+ v)))
#   }
# 
#   return(sim.phi1)
# }

# condi.phi2 <- function(y, theta, phi1, v) {
#   y0 <- 6
#   x <- y0
#   delta <- rep(NA, M) 
#   for (i in 1: M) {
#     if (x > -theta) {
#       delta[i] <- 1
#       x <- y[i]
#     } else {
#       delta[i] <- 2
#       x <- y[i]
#     }
#   }
#   idx2 <- which(delta %in% 2)
#   n2 <- length(idx2)
#   previdx2 <- idx2 - rep(1, n2)
#   ymodel2 <- rep(NA, n2) 
#   for (i in 1: n2) {  
#     ymodel2[i] <- y[idx2[i]]
#   }
#   C <- 0
#   if (y0 > -theta) {
#     C <- ymodel2[1]*y0
#     previdx2 <- previdx2[2: n2]
#     C1 <- (C + sum(ymodel2[2:n2]*y[previdx2]))/(sum(y[previdx2]^2)+y0+v)
#     sim.phi2 <- rnorm(1, mean = C1, sd = sqrt((sum(y[previdx2]^2)+y0 + v)))
#   } else {
#     C1 <- sum(ymodel1*y[previdx2])/(sum(y[previdx2]^2)+v)
#     sim.phi2 <- rnorm(1, mean = C1, sd = sqrt((sum(y[previdx2]^2)+ v)))
#   }
#   return(sim.phi2)
# }



# # saving all iterations
# Theta <- theta0
# Phi1 <- phi10
# Phi2 <- phi20
# V <- v0

#   # Gibbs step always accepted
#   vnew <- condi.v(y, theta0, phi10, phi20)
#   phi1new <- condi.phi1(y, theta0, phi20, vnew)
#   
#   phi2new <- condi.phi2(y, theta0, phi1new, vnew)
#  
#   thnew <- theta0 + rnorm(1, mean = 0, sd = 1)  # sampling new proposal value



# noraml can be outside the range. If outside the range, reject it
#   
#   if (thnew < -a || thnew > a) {
#     theta0 <- theta0
#   }  
#   else {
#     acceptance <- min(c(1, exp((-1/(2*vnew))*(sumsquare(y, thnew, phi1new, phi2new)-sumsquare(y, theta0, phi1new, phi2new)))))
#     u <- runif(1)
#     if (u < acceptance) {
#       theta0 <- thnew
#       accept <- accept + 1
#     }
#   }

#   # update candidate generating density variances
#   if(i<100){
#     new.phi1.sig<-.01
#     new.phi2.sig<-.01
#     new.theta.sig<-.01
#     new.v.sig<-.1
#   } else{
#     new.phi1.sig<-sd(phi1[1:(i-1)])+.001
#     new.phi2.sig<-sd(phi2[1:(i-1)])+.001
#     new.theta.sig<-sd(theta[1:(i-1)])+.001
#     new.v.sig<-sd(v[1:(i-1)])+.001
#   }