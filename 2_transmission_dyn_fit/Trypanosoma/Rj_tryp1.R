library('fitR')
library('MASS')
library(tidyr)
library(ggplot2)
library('coda')
library(lattice)
source('Fuctions.R')
source('tryp_det.R')
load("Tryps.RData")

Tryp_Haddon <- fix_Tryp_Haddon
Tryp_Haddon<-Tryp_Haddon[complete.cases(Tryp_Haddon),]

init.state <- c(Sv=500,Iv=1,S=80,I=1) 
theta <- c(beta1=5e-1,beta2=1e-3,gamma=0.25)

data1 <- Tryp_Haddon

posTdC <- function(theta){
  
  my_fitmodel1 <- tryp_det
  my_init.state <- init.state
  return(dLogPos6(fitmodel = my_fitmodel1,
                       theta = theta,
                       init.state =  my_init.state,
                       data = data1))
  
}

init.theta <- theta
lower <- c(beta1=1e-1,beta2=0,gamma=0)
upper <- c(beta1=1,beta2=1e-1,gamma=1)


n.iterations <- 1000
adapt.size.start <- 100
adapt.size.cooling <- 0.999
adapt.shape.start <- 200

my_mcmc.TdC <- mcmcMH(target = posTdC,
                      init.theta = theta,
                      limits = list(lower = lower,upper = upper),
                      n.iterations = n.iterations,
                      adapt.size.start = adapt.size.start,
                      adapt.size.cooling = adapt.size.cooling,
                      adapt.shape.start = adapt.shape.start)


theta <- colMeans(my_mcmc.TdC$trace[100:1000,1:3])
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

save(my_mcmc.TdC,my_mcmc.TdC2,file='Tryp.ch1.RData')