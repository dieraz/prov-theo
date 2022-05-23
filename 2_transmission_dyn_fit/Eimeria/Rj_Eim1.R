library('fitR')
library('MASS')
library(tidyr)
library(ggplot2)
library('coda')
library(lattice)
source('Functions.R')
source('Eim_det.R')
load("Eimeria.RData")

prevE_HaddonW <- fix_prevE_Haddon
prevE_HaddonW<-prevE_HaddonW[complete.cases(prevE_HaddonW), ]
prevE_HaddonW$obs <- prevE_HaddonW$obs*prevE_HaddonW$num

init.state <- c(S=80,I=1,E=0) 
theta <- c(delta=2e-6,gamma0=2e-1)

data1 <- prevE_HaddonW

posTdC <- function(theta){
  
  my_fitmodel1 <- Eim_det
  my_init.state <- init.state
  return(dLogPos6(fitmodel = my_fitmodel1,
                  theta = theta,
                  init.state =  my_init.state,
                  data = data1))
  
}

set.seed(54)
init.theta <- theta
lower <- c(delta=1e-7,gamma0=1e-4)
upper <- c(delta=1e-4,gamma0=5e-1)

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


theta <- colMeans(my_mcmc.TdC$trace[100:1000,1:2])
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

save(my_mcmc.TdC,my_mcmc.TdC2,file='Eim.ch1.RData')
