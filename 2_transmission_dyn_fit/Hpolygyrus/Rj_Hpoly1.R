library('fitR')
library('MASS')
library(tidyr)
library(ggplot2)
library('coda')
library(lattice)
source('Functions.R')
source('Hpoly_det.R')
load("Hpoly.RData")

Hpoly_HaddonW <- fix_hpoly_Haddon
Hpoly_HaddonW<-Hpoly_HaddonW[complete.cases(Hpoly_HaddonW), ]


init.state <- c(H=80,L=0,P=1) 
theta <- c(beta=4e-5,vh=9e-2)

data1 <- Hpoly_HaddonW

posTdC <- function(theta){
  
  my_fitmodel1 <- Hpoly_det
  my_init.state <- init.state
  return(dLogPos4(fitmodel = my_fitmodel1,
                  theta = theta,
                  init.state =  my_init.state,
                  data = data1))
  
}

init.theta <- theta
lower <- c(beta=1e-6,vh=1e-3)
upper <- c(beta=1e-1,vh=1)

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

save(my_mcmc.TdC,my_mcmc.TdC2,file='Hpoly.ch1.RData')