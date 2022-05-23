simulate <- function(theta,init.state,times) {
  ode <-function(t,x,params){
    
    S<-x[1]
    I<-x[2]
    E<-x[3]

    
    with(as.list(params),{
      N <- S+I  # total host population
      b <- 1.24
      K <- 42.12
      mu <- 0.0727
      alpha <- 1
      ve <- 0
      taue <- 717
      phie <- 1/9
      
      
      dS <- b*N*(K-N)/K + gamma0*I - alpha*delta*S*E - mu*S
      dI <- alpha*delta*S*E - (mu + ve+gamma0)*I
      dE <- taue*I - phie*E


      dx <- c(dS,dI,dE) 
      list(dx)
    })
  }
  
  
  traj <- as.data.frame(lsoda(init.state, times, ode, theta))
  
  return(traj)
  }

rPointObs <- function(model.point, theta){
  
  obs.point <- rpois(n=1, lambda=model.point[["I"]])  
  return(c(obs=obs.point))
}

dPointObs <- function(data.point, model.point, theta,log = FALSE){
  
  return(dpois(x=data.point[["obs"]],lambda=model.point[["I"]],log=log))
  
}

dprior <- function(theta, log = FALSE) {
  
  log.prior.delta <- dunif(theta[["delta"]], min = 1e-7, max = 1e-4, log = TRUE)
  log.prior.gamma0 <- dpois(round(1/theta[["gamma0"]]),lambda = 2.8, log = TRUE)
  log.sum = log.prior.delta + log.prior.gamma0
  
  
  return(ifelse(log, log.sum, exp(log.sum)))
  
}

name <- "Coinfection model"
state.names <- c("S","I","E")
theta.names <- c("delta","gamma0")

Eim_det <- fitmodel(name, state.names, theta.names,simulate, rPointObs, dprior,dPointObs)  
