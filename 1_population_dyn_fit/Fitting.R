library('fitR')
library('MASS')
library(tidyr)
library(ggplot2)
library('coda')
library(lattice)

source('Functions.R')
source('popdyn_det.R')

load("POP_Haddon.RData")

POP_Haddon <- fix_POP_Haddon
POP_Haddon<-POP_Haddon[complete.cases(POP_Haddon), ]


init.state <- c(N=50) 
theta <- c(b=0.5,K=50,mu=1/24)

data1 <- POP_Haddon

posTdC <- function(theta){
  
  my_fitmodel1 <- Popdyn_det
  my_init.state <- init.state
  return(dLogPosPOPDYN(fitmodel1 = my_fitmodel1,
                       theta = theta,
                       init.state =  my_init.state,
                       data1 = data1))
  
}


init.theta <- theta
lower <- c(b=0.5,K=25,mu=1/52)
upper <- c(b=20,K=100,mu=1/8)


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

theta <- colMeans(my_mcmc.TdC$trace[100:1000,1:3])
covmat <- my_mcmc.TdC$covmat.empirical

n.iterations <- 100000
adapt.size.start <- 200
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

save(my_mcmc.TdC2,file='POPDYN.ch1.RData')