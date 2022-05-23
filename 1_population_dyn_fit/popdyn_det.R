simulate <- function(theta,init.state,times) {
  ode <-function(t,x,params){
    
    N <- x[1]
    
    
    with(as.list(params),{

      dN <- b*N*(K-N)/K  - mu*N

      dx <- c(dN) 
      list(dx)
    })
  }
  
  traj <- as.data.frame(lsoda(init.state, times, ode, theta))
  
  return(traj)
}



rPointObs <- function(model.point, theta){
  
  obs.point <- rpois(n=1, lambda=model.point[["N"]])  
  return(c(obs=obs.point))
}

dPointObs <- function(data.point, model.point, theta,log = FALSE){
  
  return(dpois(x=data.point[["obs"]],lambda=model.point[["N"]],log=log))
  
}

dprior <- function(theta, log = FALSE) {
  log.prior.b <- dunif(theta[["b"]], min = 0.5, max = 2, log = TRUE)  
  log.prior.K <- dunif(theta[["K"]], min = 25, max = 100, log = TRUE)  
  log.prior.mu <- dunif(theta[["mu"]], min = 1/52, max = 1/8, log = TRUE)  

  log.sum = log.prior.b + log.prior.K + log.prior.mu
  
  
  return(ifelse(log, log.sum, exp(log.sum)))
  
}

name <- "Population dynamics model"
state.names <- c("N")
theta.names <- c("b","K","mu")

Popdyn_det <- fitmodel(name, state.names, theta.names,simulate, rPointObs, dprior,dPointObs)  
