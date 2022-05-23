simulate <- function(theta,init.state,times) {
  ode <-function(t,x,params){
    
    H<-x[1]
    L<-x[2]
    P<-x[3]
    Ldet <- x[4]

    
    with(as.list(params),{
      b <- 1.24
      K <- 42.12
      mu <- 0.0727
      alpha <- 1
      epsilon <- 0
      k <- 0.1
      tauh <- 42
      phih <- 1/8

      dH <- (b*(K-H)/K - mu)*H - vh*P
      dP <- beta*L*H - (mu + vh + epsilon)*P - (P^2*vh*(k+1))/(H*k)
      dL <- tauh*P - phih*L
      dLdet <- tauh*P

      dx <- c(dH,dL,dP,dLdet) 
    list(dx)
    })
  }
  
  # put incidence at 0 in init.state
  init.state["Ldet"] <- 0
  
  traj <- as.data.frame(lsoda(init.state, times, ode, theta))
  
  traj$Ldet <- c(0,diff(traj$Ldet))
  
  return(traj)
  }

rPointObs <- function(model.point, theta){
  ##only model
  obs.point <- rnbinom(n=1,size=1,mu=model.point[["Ldet"]])
  return(c(obs=obs.point))
}

dPointObs <- function(data.point, model.point, theta, log = FALSE){
  ##model against data
  return(dnbinom(x=data.point[["obs"]],size=1,mu=model.point[["Ldet"]],log=log))
}

dprior <- function(theta, log = FALSE) {
  
  log.prior.beta <- dunif(theta[["beta"]], min = 1e-6, max = 1e-1, log = TRUE)
  log.prior.vh <- dunif(theta[["vh"]], min = 0, max = 1, log = TRUE)

  
  log.sum = log.prior.beta + log.prior.vh
  
  
  return(ifelse(log, log.sum, exp(log.sum)))
  
}

name <- "Coinfection model"
state.names <- c("H","L","P")
theta.names <- c("beta","vh")

Hpoly_det <- fitmodel(name, state.names, theta.names,simulate, rPointObs, dprior,dPointObs)  
