simulate <- function(theta,init.state,times) {
  ode <-function(t,x,params){
    
    Sm <- x[1]
    Am <- x[2]
    Lm <-x[3]
    Sf <- x[4]
    Af <-x[5]
    Lf <-x[6]
    
    
    with(as.list(params),{
      N <- Sm+Am+Lm+Sf+Af+Lf  # total host population
      b <- 1.24
      K <- 42.12
      mu <- 0.0727
      nu <- 0
      
      
      dSm <- b*N/2*(K-N)/K - beta1*Sm*Am - beta2*Sm*Af - mu*Sm
      dAm <- beta1*Sm*Am + beta2*Sm*Af + eta*Lm - (epsilon+mu+nu)*Am
      dLm <- epsilon*Am - (eta+mu)*Lm
      
      dSf <- b*N/2*(K-N)/K - beta3*Sf*Am - beta4*Sf*Af - mu*Sf
      dAf <- beta3*Sf*Am + beta4*Sf*Af + eta*Lf - (epsilon+mu+nu)*Af 
      dLf <- epsilon*Af - (eta+mu)*Lf
      
      
      dx <- c(dSm,dAm,dLm,dSf,dAf,dLf) 
      list(dx)
    })
  }
  
  traj <- as.data.frame(lsoda(init.state, times, ode, theta))
  traj$prev <-(traj$Am +traj$Lm)/(traj$Sm + traj$Am +traj$Lm)
  
  return(traj)
}

rPointObs <- function(model.point, theta){
  
  obs.point <- rbinom(1,round(model.point[["Sm"]]+model.point[["Am"]]+model.point[["Lm"]]),model.point[["prev"]])/round(model.point[["Sm"]]+model.point[["Am"]]+model.point[["Lm"]])
  return(c(obs=obs.point))
}

dPointObs <- function(data.point, model.point, theta,log = FALSE){
  size <-  round(data.point[["num"]])
  x <- round(data.point[["obs"]]*data.point[["num"]])
  prob <- model.point[["prev"]]
  return(dbinom(x=x,size=size,prob = prob,log=log))
  
}


dprior <- function(theta, log = FALSE) {
  
  log.prior.beta1 <- dunif(theta[["beta1"]], min = 1e-7, max = 1e-1, log = TRUE)  
  log.prior.beta2 <- dunif(theta[["beta2"]], min = 1e-7, max = 1e-1, log = TRUE)  
  log.prior.beta3 <- dunif(theta[["beta3"]], min = 1e-7, max = 1e-1, log = TRUE)  
  log.prior.beta4 <- dunif(theta[["beta4"]], min = 1e-7, max = 1e-1, log = TRUE)
  log.prior.eta <- dunif(theta[["eta"]], min = 1e-7, max = 1, log = TRUE)
  log.prior.epsilon <- dunif(theta[["epsilon"]], min = 1e-7, max = 1, log = TRUE)
  
  log.sum = log.prior.beta1 + log.prior.beta2 + log.prior.beta3 + log.prior.beta4 + log.prior.eta + log.prior.epsilon
  
  
  return(ifelse(log, log.sum, exp(log.sum)))
  
}

name <- "Herpes model"
state.names <- c("Sm","Am","Lm","Sf","Af","Lf")
theta.names <- c("beta1","beta2","beta3","beta4","eta","epsilon")

herpes_det_m <- fitmodel(name, state.names, theta.names,simulate, rPointObs, dprior,dPointObs)  
