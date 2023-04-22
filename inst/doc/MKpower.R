## -----------------------------------------------------------------------------
library(MKpower)

## -----------------------------------------------------------------------------
## with continuity correction
power.prop1.test(p1 = 0.4, power = 0.8)
## without continuity correction
power.prop1.test(p1 = 0.4, power = 0.8, cont.corr = FALSE)

## -----------------------------------------------------------------------------
## with continuity correction
ssize.propCI(prop = 0.4, width = 0.14)
## without continuity correction
ssize.propCI(prop = 0.4, width = 0.14, method = "wald")

## -----------------------------------------------------------------------------
## identical results as power.t.test, since sd = sd1 = sd2 = 1
power.welch.t.test(n = 20, delta = 1)
power.welch.t.test(power = .90, delta = 1)
power.welch.t.test(power = .90, delta = 1, alternative = "one.sided")

## sd1 = 0.5, sd2 = 1
power.welch.t.test(delta = 1, sd1 = 0.5, sd2 = 1, power = 0.9)

## -----------------------------------------------------------------------------
## slightly more conservative than Welch t-test
power.hsu.t.test(n = 20, delta = 1)
power.hsu.t.test(power = .90, delta = 1)
power.hsu.t.test(power = .90, delta = 1, alternative = "one.sided")

## sd1 = 0.5, sd2 = 1
power.welch.t.test(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9)
power.hsu.t.test(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9)

## -----------------------------------------------------------------------------
## 3 groups
cbind(rep(1,2), -diag(2))
## 4 groups
cbind(rep(1,3), -diag(3))

## -----------------------------------------------------------------------------
## Example from Maxwell and Delaney (2004) according to Shieh (2020)
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.8)
power.ancova(n = rep(45/3, 3), mu=c(7.5366, 11.9849, 13.9785), var = 29.0898)
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.9)
power.ancova(n = rep(57/3, 3), mu=c(7.5366, 11.9849, 13.9785), var = 29.0898)

## -----------------------------------------------------------------------------
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.8,
             nr.covs = 2)
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.8,
             nr.covs = 3)
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.8,
             nr.covs = 5)
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.8,
             nr.covs = 10)

## -----------------------------------------------------------------------------
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.8, 
             group.ratio = c(1, 1.25, 1.5))
power.ancova(n = c(13, 16, 19), mu=c(7.5366, 11.9849, 13.9785), var = 29.0898,  
             group.ratio = c(1, 1.25, 1.5))
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.8, 
             group.ratio = c(1, 0.8, 2/3))
power.ancova(n = c(17, 14, 12), mu=c(7.5366, 11.9849, 13.9785), var = 29.0898,  
             group.ratio = c(1, 0.8, 2/3))

## -----------------------------------------------------------------------------
rx <- function(n) rnorm(n, mean = 0, sd = 1) 
ry <- function(n) rnorm(n, mean = 0.5, sd = 1) 
## two-sample
sim.ssize.wilcox.test(rx = rx, ry = ry, n.max = 100, iter = 1000)
sim.ssize.wilcox.test(rx = rx, ry = ry, n.min = 65, n.max = 70, 
                      step.size = 1, iter = 1000, BREAK = FALSE)
power.t.test(delta = 0.5, power = 0.8)

## one-sample
sim.ssize.wilcox.test(rx = ry, n.max = 100, iter = 1000, type = "one.sample")
sim.ssize.wilcox.test(rx = ry, type = "one.sample", n.min = 33, 
                      n.max = 38, step.size = 1, iter = 1000, BREAK = FALSE)
power.t.test(delta = 0.5, power = 0.8, type = "one.sample")

## two-sample
rx <- function(n) rgamma(n, scale = 2, shape = 1.5)
ry <- function(n) rgamma(n, scale = 4, shape = 1.5) # equal shape
## two-sample
sim.ssize.wilcox.test(rx = rx, ry = ry, n.max = 100, iter = 1000)
sim.ssize.wilcox.test(rx = rx, ry = ry, n.min = 25, n.max = 30, 
                      step.size = 1, iter = 1000, BREAK = FALSE)

## -----------------------------------------------------------------------------
## examples from Table III in Zhu and Lakkis (2014)
power.nb.test(mu0 = 5.0, RR = 2.0, theta = 1/0.5, duration = 1, power = 0.8, approach = 1)
power.nb.test(mu0 = 5.0, RR = 2.0, theta = 1/0.5, duration = 1, power = 0.8, approach = 2)
power.nb.test(mu0 = 5.0, RR = 2.0, theta = 1/0.5, duration = 1, power = 0.8, approach = 3)

## -----------------------------------------------------------------------------
## see n2 on page 1202 of Chu and Cole (2007)
power.diagnostic.test(sens = 0.99, delta = 0.14, power = 0.95) # 40
power.diagnostic.test(sens = 0.99, delta = 0.13, power = 0.95) # 43
power.diagnostic.test(sens = 0.99, delta = 0.12, power = 0.95) # 47

## -----------------------------------------------------------------------------
## see Table 2 of Dobbin et al. (2008)
g <- 0.1
fc <- 1.6
ssize.pcc(gamma = g, stdFC = fc, nrFeatures = 22000)

## -----------------------------------------------------------------------------
## functions to simulate the data
## group 1
rx <- function(n) rnorm(n, mean = 0, sd = 1)
rx.H0 <- rx
## group 2
ry <- function(n) rnorm(n, mean = 3, sd = 3)
ry.H0 <- function(n) rnorm(n, mean = 0, sd = 3)
## theoretical results
power.welch.t.test(n = 10, delta = 3, sd1 = 1, sd2 = 3)
power.hsu.t.test(n = 10, delta = 3, sd1 = 1, sd2 = 3)
## simulation
res.t <- sim.power.t.test(nx = 10, rx = rx, rx.H0 = rx.H0,
                          ny = 10, ry = ry, ry.H0 = ry.H0)
res.t
res.w <- sim.power.wilcox.test(nx = 10, rx = rx, rx.H0 = rx.H0,
                               ny = 10, ry = ry, ry.H0 = ry.H0)
res.w

## ---- fig.width=7, fig.height=7-----------------------------------------------
## t-tests
hist(res.t)
qqunif(res.t, alpha = 0.05)
volcano(res.t, hex = TRUE)
##  Wilcoxon-Mann-Whitney tests
hist(res.w)
qqunif(res.w, alpha = 0.05)

## -----------------------------------------------------------------------------
sessionInfo()

