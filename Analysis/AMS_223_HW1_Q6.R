# AMS 223 Time Series HW1 Q6
# Cheng-Han Yu, Dept of Statistics UC Santa Cruz
# Time series and Baysian inference overview

## ---- Q6a1 ----
# (a) plot the data
# EEG at channel F3
eegf3 <- read.table("http://users.soe.ucsc.edu/~raquel/tsbook/data/eegF3.dat")
eegf3 <- as.vector(t(eegf3))
# head(eegf3)  # check data
plot(eegf3, type = "l", axes = F, xlab = "time", ylab = " ")
axis(1)
axis(2)

## ---- Q6a2 ----
# GDP
gdp <- read.table("http://users.soe.ucsc.edu/~raquel/tsbook/data/gdp.dat",
                  header = TRUE, skip = 2)
# head(gdp)  # check data
par(mfrow = c(3, 3))
for (i in 2: ncol(gdp)) {
  plot(gdp[, 1], gdp[, i], type = "l", xlab = "year", ylab = " ", axes = F,
       main = colnames(gdp)[i], cex.axis = 0.5)
  axis(1)
  axis(2)
}
dev.off()

## ---- Q6a3 ----
# Southern Oscillation Index (SOI)
soi <- read.table("http://users.soe.ucsc.edu/~raquel/tsbook/data/soi.dat")
# head(soi)  # check data
ts.plot(soi, gpars = list(xlab = "time", ylab = " ", axes = F))
axis(1)
axis(2)

## ---- Q6b11 ----
# (b) ACF, smoothing SOI, differencing GDP
# (1) ACF
acf(eegf3)

## ---- Q6b12 ----
par(mfrow = c(3,3))
for (i in 2: ncol(gdp)) {
  acf(gdp[, i], main = colnames(gdp)[i])
}

## ---- Q6b13 ----
acf(soi)

## ---- Q6b21 ----
# (2) smoothing SOI
# moving avg of order 5 with equal weights
soi.ma5 <- filter(soi, filter = c(.2, .2, .2, .2, .2), side = 2)
# moving avg of order 11 with equal weights
soi.ma11 <- filter(soi, filter = rep(1, 11)/11, side = 2)
par(mfrow = c(2, 1))
ts.plot(soi.ma5, 
        gpars = list(xlab = "time", ylab = " ", axes = F, cex.main = 0.95),
        main = "Smoothing Series of SOI: MA order 5 with equal weights")
axis(1)
axis(2)
ts.plot(soi.ma11, 
        gpars = list(xlab = "time", ylab = " ", axes = F, cex.main = 0.95),
        main = "Smoothing Series of SOI: MA order 11 with equal weights")
axis(1)
axis(2)


## ---- Q6b22 ----
# moving avg of order 5 with unequal weights (0.1, 0.15, 0.2, 0.25, 0.3)
soi.ma5un1 <- filter(soi, filter = c(0.1, 0.15, 0.2, 0.25, 0.3), side = 2)
# moving avg of order 5 with unequal weights (0.3, 0.25, 0.2, 0.15, 0.1)
soi.ma5un2 <- filter(soi, filter = c(0.3, 0.25, 0.2, 0.15, 0.1), side = 2)
par(mfrow = c(2, 1))
ts.plot(soi.ma5un1, gpars = list(xlab = "time", ylab = " ", axes = F, 
                                 cex.main = 0.9),
        main = "Smoothing SOI: MA order 5 with unequal weights 
        (0.1, 0.15, 0.2, 0.25, 0.3)")
axis(1)
axis(2)
ts.plot(soi.ma5un2, gpars = list(xlab = "time", ylab = " ", axes = F,
                                 cex.main = 0.9),
        main = "Smoothing SOI: MA order 5 with unequal weights
        (0.3, 0.25, 0.2, 0.15, 0.1)")
axis(1)
axis(2)

## ---- Q6b31 ----
# (3) differencing GDP ACF
gdp.d1 <- matrix(NA, nrow = nrow(gdp) - 1, ncol = ncol(gdp) - 1)
for (i in 2: ncol(gdp)) {
  gdp.d1[, i - 1] <- diff(gdp[, i])
}
gdp.d1 <- data.frame(cbind(gdp[(2:nrow(gdp)), 1], gdp.d1))
colnames(gdp.d1) = paste("1st dif", colnames(gdp))

par(mfrow = c(3, 3))
for (i in 2: ncol(gdp)) {
  plot(gdp.d1[, 1], gdp.d1[, i], type = "l", xlab = "year", ylab = " ", axes = F,
       main = colnames(gdp)[i])
  axis(1)
  axis(2)
}


## ---- Q6b32 ----
par(mfrow = c(3,3))
for (i in 2: ncol(gdp)) {
  acf(gdp.d1[, i], main = paste("ACF of 1st diff of", colnames(gdp)[i]),
      cex.main = 0.85)
}

## ---- Q6b33 ----
gdp.d2 <- matrix(NA, nrow = nrow(gdp) - 2, ncol = ncol(gdp) - 1)
for (i in 2: ncol(gdp)) {
  gdp.d2[, i - 1] <- diff(gdp[, i], difference = 2)
}
gdp.d2 <- data.frame(cbind(gdp[(3: nrow(gdp)), 1], gdp.d2))
colnames(gdp.d2) = paste("2nd dif", colnames(gdp))
par(mfrow = c(3, 3))
for (i in 2: ncol(gdp)) {
  plot(gdp.d2[, 1], gdp.d2[, i], type = "l", xlab = "year", ylab = " ", axes = F,
       main = colnames(gdp)[i])
  axis(1)
  axis(2)
}


## ---- Q6b34 ----
par(mfrow = c(3,3))
for (i in 2: ncol(gdp)) {
  acf(gdp.d2[, i], main = paste("ACF of 2nd diff of", colnames(gdp)[i]), 
      cex.main = 0.85)
}


