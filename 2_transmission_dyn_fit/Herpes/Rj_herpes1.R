library('fitR')
library('MASS')
library(tidyr)
library(ggplot2)
library('coda')
library(lattice)
source('Functions.R')
source('herpes_det_ffixed.R')
source('herpes_det_mfixed.R')
load("Herpes.RData")

WMHV_Haddon_F <- fix_HV_Haddon_F
WMHV_Haddon_F<-WMHV_Haddon_F[complete.cases(WMHV_Haddon_F), ]

WMHV_Haddon_M <- fix_HV_Haddon_M
WMHV_Haddon_M<-WMHV_Haddon_M[complete.cases(WMHV_Haddon_M), ]

init.state <- c(Sm=25,Am=1,Lm=0,Sf=25,Af=0,Lf=0) 
theta <- c(beta1=1e-2,beta2=6e-3,beta3=2e-3,beta4=2e-3,eta=0.25,epsilon=0.25)

data1 <- WMHV_Haddon_F
data2 <- WMHV_Haddon_M


posTdC <- function(theta){
  
  my_fitmodel1 <- herpes_det_f
  my_fitmodel2<- herpes_det_m
  my_init.state <- init.state
  return(dLogPosherpes(fitmodel1 = my_fitmodel1,
                       fitmodel2 = my_fitmodel2,
                       theta = theta,
                       init.state =  my_init.state,
                       data1 = data1,
                       data2 = data2))
  
}

init.theta <- theta
lower <- c(beta1=0,beta2=0,beta3=0,beta4=0,eta=0,epsilon=0)
upper <- c(beta1=1e-1,beta2=1e-1,beta3=1e-1,beta4=1e-1,eta=1,epsilon=1)


n.iterations <- 1000
adapt.size.start <- 100
adapt.size.cooling <- 0.999
adapt.shape.start <- 200

my_mcmc.TdC <- mcmcMH(target = posTdC,
                      init.theta = theta,
                      #proposal.sd = proposal.sd,
                      limits = list(lower = lower,upper = upper),
                      n.iterations = n.iterations,
                      adapt.size.start = adapt.size.start,
                      adapt.size.cooling = adapt.size.cooling,
                      adapt.shape.start = adapt.shape.start)

theta <- colMeans(my_mcmc.TdC$trace[100:1000,1:6])
covmat <- my_mcmc.TdC$covmat.empirical

n.iterations <- 100000
adapt.size.start <- 100
adapt.size.cooling <- 0.999
adapt.shape.start <- 200

my_mcmc.TdC2 <- mcmcMH(target = posTdC,
                       init.theta = theta,
                       covmat = covmat,
                       limits = list(lower = lower,upper = upper),
                       n.iterations = n.iterations,
                       adapt.size.start = adapt.size.start,
                       adapt.size.cooling = adapt.size.cooling,
                       adapt.shape.start = adapt.shape.start)

save(my_mcmc.TdC,my_mcmc.TdC2,file='WMHV.ch1.RData')