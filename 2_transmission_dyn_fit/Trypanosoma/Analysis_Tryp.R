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

load("Tryps.RData")
source('tryp_det.R')
source('Functions.R')

load("Tryp.ch1.RData")
my_mcmc.TdC2$acceptance.rate

trace1 <- mcmc(my_mcmc.TdC2$trace)
trace.burn1 <- burnAndThin(trace1, burn = 5000)
trace.burn.thin1 <- burnAndThin(trace.burn1, thin = 10)
theta1 <- colMeans(my_mcmc.TdC2$trace[,])

load("Tryp.ch2.RData")
my_mcmc.TdC2$acceptance.rate

trace2 <- mcmc(my_mcmc.TdC2$trace)
trace.burn2 <- burnAndThin(trace2, burn = 5000)
trace.burn.thin2 <- burnAndThin(trace.burn2, thin = 10)
theta2 <- colMeans(my_mcmc.TdC2$trace[,])

load("Tryp.ch3.RData")
my_mcmc.TdC2$acceptance.rate

trace3 <- mcmc(my_mcmc.TdC2$trace)
trace.burn3 <- burnAndThin(trace3, burn = 5000)
trace.burn.thin3 <- burnAndThin(trace.burn3, thin = 10)
theta3 <- colMeans(my_mcmc.TdC2$trace[,])

trace.info <- mcmc.list(trace.burn.thin1,trace.burn.thin2,trace.burn.thin3)
densityplot(trace.info)

gelman.diag(trace.info)

trace.combined <- ldply(trace.info)
theta.bar <- colMeans(trace.combined[tryp_det$theta.names])

Tryp_HaddonW <- fix_Tryp_Haddon
Tryp_HaddonW<-Tryp_HaddonW[complete.cases(Tryp_HaddonW), ]

init.state <- c(Sv=500,Iv=1,S=80,I=1) 
log.like.theta.bar1 <- dTrajObs_DE6(tryp_det, theta.bar, init.state, data = Tryp_HaddonW, 
                                     log = TRUE)
D.theta.bar1 <- -2 * log.like.theta.bar1
p.D <- var(-2 * trace.combined$log.density)/2
DIC1 <- D.theta.bar1 + 2 * p.D

gelman.diag(trace.info)

plotFit_DE(tryp_det,theta.bar,init.state,data = Tryp_HaddonW,n.replicates = 20)
