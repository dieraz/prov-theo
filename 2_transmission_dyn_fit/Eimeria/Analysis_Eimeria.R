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

load("Eimeria.RData")
source('Eim_det.R')
source('Functions.R')

load("Eim.ch1.RData")
my_mcmc.TdC2$acceptance.rate

trace1 <- mcmc(my_mcmc.TdC2$trace)
trace.burn1 <- burnAndThin(trace1, burn = 5000)
trace.burn.thin1 <- burnAndThin(trace.burn1, thin = 10)
theta1 <- colMeans(my_mcmc.TdC2$trace[,])

load("Eim.ch2.RData")
my_mcmc.TdC2$acceptance.rate

trace2 <- mcmc(my_mcmc.TdC2$trace)
trace.burn2 <- burnAndThin(trace2, burn = 5000)
trace.burn.thin2 <- burnAndThin(trace.burn2, thin = 10)
theta2 <- colMeans(my_mcmc.TdC2$trace[,])

load("Eim.ch3.RData")
my_mcmc.TdC2$acceptance.rate

trace3 <- mcmc(my_mcmc.TdC2$trace)
trace.burn3 <- burnAndThin(trace3, burn = 5000)
trace.burn.thin3 <- burnAndThin(trace.burn3, thin = 10)
theta3 <- colMeans(my_mcmc.TdC2$trace[,])

trace.info <- mcmc.list(trace.burn.thin1,trace.burn.thin2,trace.burn.thin3)
densityplot(trace.info)

gelman.diag(trace.info)

trace.combined <- ldply(trace.info)
theta.bar <- colMeans(trace.combined[Eim_det$theta.names])

prevE_HaddonW <- fix_prevE_Haddon
prevE_HaddonW<-prevE_HaddonW[complete.cases(prevE_HaddonW), ]
prevE_HaddonW$obs <- prevE_HaddonW$obs*prevE_HaddonW$num

init.state <- c(S=50,I=1,R=0)

log.like.theta.bar1 <- dTrajObs_DE6(Eim_det, theta.bar, init.state, data = prevE_HaddonW, 
                                     log = TRUE)
D.theta.bar1 <- -2 * log.like.theta.bar1
p.D <- var(-2 * trace.combined$log.density)/2
DIC1 <- D.theta.bar1 + 2 * p.D

gelman.diag(trace.info)

plotFit_DE(Eim_det,theta.bar,init.state,data = prevE_HaddonW,n.replicates = 20)
