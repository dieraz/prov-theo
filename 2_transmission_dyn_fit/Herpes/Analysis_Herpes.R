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

load("herpes.RData")
source('Functions.R')
source('herpes_det_ffixed.R')
source('herpes_det_mfixed.R')
init.state <- c(Sm=25,Am=1,Lm=0,Sf=25,Af=0,Lf=0) 

load("WMHV.ch1.RData")
my_mcmc.TdC2$acceptance.rate

trace1 <- mcmc(my_mcmc.TdC2$trace)
trace.burn1 <- burnAndThin(trace1, burn = 5000)
trace.burn.thin1 <- burnAndThin(trace.burn1, thin = 10)
theta1 <- colMeans(my_mcmc.TdC2$trace[,])

load("WMHV.ch2.RData")
my_mcmc.TdC2$acceptance.rate

trace2 <- mcmc(my_mcmc.TdC2$trace)
trace.burn2 <- burnAndThin(trace2, burn = 5000)
trace.burn.thin2 <- burnAndThin(trace.burn2, thin = 10)
theta2 <- colMeans(my_mcmc.TdC2$trace[,])

load("WMHV.ch3.RData")
my_mcmc.TdC2$acceptance.rate

trace3 <- mcmc(my_mcmc.TdC2$trace)
trace.burn3 <- burnAndThin(trace3, burn = 5000)
trace.burn.thin3 <- burnAndThin(trace.burn3, thin = 10)
theta3 <- colMeans(my_mcmc.TdC2$trace[,])


trace.info <- mcmc.list(trace.burn.thin1,trace.burn.thin2,trace.burn.thin3)
densityplot(trace.info)

gelman.diag(trace.info)

WMHV_Haddon_F <- fix_HV_Haddon_F
WMHV_Haddon_F<-WMHV_Haddon_F[complete.cases(WMHV_Haddon_F), ]

WMHV_Haddon_M <- fix_HV_Haddon_M
WMHV_Haddon_M<-WMHV_Haddon_M[complete.cases(WMHV_Haddon_M), ]

trace.combined <- ldply(trace.info)
theta.bar <- colMeans(trace.combined[herpes_det_f$theta.names])

log.like.theta.bar1 <- dTrajObs_DE_hf(herpes_det_f, theta.bar, init.state, data = WMHV_Haddon_F, 
                                     log = TRUE)
D.theta.bar1 <- -2 * log.like.theta.bar1
p.D <- var(-2 * trace.combined$log.density)/2
DIC1 <- D.theta.bar1 + 2 * p.D

log.like.theta.bar2 <- dTrajObs_DE_hm(herpes_det_m, theta.bar, init.state, data = WMHV_Haddon_M, 
                                     log = TRUE)
D.theta.bar2 <- -2 * log.like.theta.bar2
DIC2 <- D.theta.bar2 + 2 * p.D

DIC <- -2*(log.like.theta.bar1+log.like.theta.bar2) + 2 * p.D

gelman.diag(trace.info)

plotFit_DE(herpes_det_f,theta.bar,init.state,data = WMHV_Haddon_F,n.replicates = 20)
plotFit_DE(herpes_det_m,theta.bar,init.state,data = WMHV_Haddon_M,n.replicates = 20)
