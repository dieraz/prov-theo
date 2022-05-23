simulate <- function(theta,init.state,times) {
  ode <-function(t,x,params){
    
    Sv <- x[1]
    Iv <- x[2]
    S <-x[3]
    I <- x[4]
    
    with(as.list(params),{
      Nv <- Sv+Iv # total host population
      N <- S+I # total host population
      b <- 1.24
      K <- 42.12
      mu <- 0.0727
      
      bv <- 100
      Kv <- K*2
      muv <- 1/4

      dSv <- bv*Nv*(Kv-Nv)/Kv - beta1*Sv*I - muv*Sv
      dIv <- beta1*Sv*I - muv*Iv
      dS <- b*N*(K-N)/K - beta2*S*Iv - mu*S + gamma*I
      dI <- beta2*S*Iv - (gamma+mu)*I
 
      dx <- c(dSv,dIv,dS,dI) 
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
  
  log.prior.beta1 <- dunif(theta[["beta1"]], min = 1e-7, max = 1, log = TRUE)  
  log.prior.beta2 <- dunif(theta[["beta2"]], min = 1e-7, max = 1, log = TRUE)  
  log.prior.gamma <- dunif(theta[["gamma"]], min = 1e-7, max = 1, log = TRUE)
  
  log.sum = log.prior.beta1 + log.prior.beta2 + log.prior.gamma
  
  
  return(ifelse(log, log.sum, exp(log.sum)))
  
}

name <- "Bart model"
state.names <- c("Sv","Iv","S","I")
theta.names <- c("beta1","beta2","gamma")

bart_det <- fitmodel(name, state.names, theta.names,simulate, rPointObs, dprior,dPointObs)  
