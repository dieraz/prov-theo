library('fitR')
library('MASS')
library(tidyr)
library(ggplot2)
library('coda')
library(lattice)
library(reshape2)
library(plyr)
library(stringr)
library(gridExtra)

load("Hpoly.RData")
source('Hpoly_det.R')
source('Functions.R')

load("Hpoly.ch1.RData")
my_mcmc.TdC2$acceptance.rate

trace1 <- mcmc(my_mcmc.TdC2$trace)
trace.burn1 <- burnAndThin(trace1, burn = 5000)
trace.burn.thin1 <- burnAndThin(trace.burn1, thin = 10)
theta1 <- colMeans(my_mcmc.TdC2$trace[,])

load("Hpoly.ch2.RData")
my_mcmc.TdC2$acceptance.rate

trace2 <- mcmc(my_mcmc.TdC2$trace)
trace.burn2 <- burnAndThin(trace2, burn = 5000)
trace.burn.thin2 <- burnAndThin(trace.burn2, thin = 10)
theta2 <- colMeans(my_mcmc.TdC2$trace[,])

load("Hpoly.ch3.RData")
my_mcmc.TdC2$acceptance.rate

trace3 <- mcmc(my_mcmc.TdC2$trace)
trace.burn3 <- burnAndThin(trace3, burn = 5000)
trace.burn.thin3 <- burnAndThin(trace.burn3, thin = 10)
theta3 <- colMeans(my_mcmc.TdC2$trace[,])

trace.info <- mcmc.list(trace.burn.thin1,trace.burn.thin2,trace.burn.thin3)
densityplot(trace.info)

gelman.diag(trace.info)

trace.combined <- ldply(trace.info)
theta.bar <- colMeans(trace.combined[Hpoly_det$theta.names])

Hpoly_HaddonW <- fix_hpoly_Haddon
Hpoly_HaddonW<-Hpoly_HaddonW[complete.cases(Hpoly_HaddonW), ]

init.state <- c(H=80,L=0,P=1)


log.like.theta.bar1 <- dTrajObs_DE4(Hpoly_det, theta.bar, init.state, data = Hpoly_HaddonW, 
                                     log = TRUE)
D.theta.bar1 <- -2 * log.like.theta.bar1
p.D <- var(-2 * trace.combined$log.density)/2
DIC1 <- D.theta.bar1 + 2 * p.D

gelman.diag(trace.info)

plotFit_DE(Hpoly_det,theta.bar,init.state,data = Hpoly_HaddonW,n.replicates = 20)


