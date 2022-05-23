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

source('popdyn_det.R')
source('Functions.R')
load("POP_Haddon.RData")
load("POPDYN.ch1.RData")

my_mcmc.TdC2$acceptance.rate

trace1 <- mcmc(my_mcmc.TdC2$trace)
trace.burn1 <- burnAndThin(trace1, burn = 5000)
trace.burn.thin1 <- burnAndThin(trace.burn1, thin = 100)
theta1 <- colMeans(trace.burn.thin1)

load("POPDYN.ch2.RData")
my_mcmc.TdC2$acceptance.rate

trace2 <- mcmc(my_mcmc.TdC2$trace)
trace.burn2 <- burnAndThin(trace2, burn = 5000)
trace.burn.thin2 <- burnAndThin(trace.burn2, thin = 100)
theta2<- colMeans(trace.burn.thin2)

load("POPDYN.ch3.RData")
my_mcmc.TdC2$acceptance.rate

trace3 <- mcmc(my_mcmc.TdC2$trace)
trace.burn3 <- burnAndThin(trace3, burn = 5000)
trace.burn.thin3 <- burnAndThin(trace.burn3, thin = 100)
theta3 <- colMeans(trace.burn.thin3)

load("POPDYN.ch4.RData")
my_mcmc.TdC2$acceptance.rate

trace4 <- mcmc(my_mcmc.TdC2$trace)
trace.burn4 <- burnAndThin(trace4, burn = 5000)
trace.burn.thin4 <- burnAndThin(trace.burn4, thin = 100)
theta4 <- colMeans(trace.burn.thin4)

trace.info <- mcmc.list(trace.burn.thin1,trace.burn.thin2,trace.burn.thin3,trace.burn.thin4)
gelman.diag(trace.info)

trace.combined <- ldply(trace.info)
theta.bar <- colMeans(trace.combined[Popdyn_det$theta.names])

POP_Haddon <- fix_POP_Haddon
POP_Haddon<-POP_Haddon[complete.cases(POP_Haddon), ]
init.state <- c(N=50) 
log.like.theta.bar <- dTrajObs_POPDYN(Popdyn_det, theta.bar, init.state, data = POP_Haddon, 
                                     log = TRUE)
D.theta.bar <- -2 * log.like.theta.bar
p.D <- var(-2 * trace.combined$log.density)/2
DIC <- D.theta.bar + 2 * p.D

plotFit_DE(Popdyn_det,theta.bar,init.state,data = POP_Haddon,n.replicates = 20)
