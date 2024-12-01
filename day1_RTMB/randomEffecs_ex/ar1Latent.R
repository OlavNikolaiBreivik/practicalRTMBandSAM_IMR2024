
library(RTMB)
y = readRDS("randomEffecs_ex/latentAR1Process.Rds")

#Set up and estimate model
data = list(y = y)
par = list(logsd = c(0,0),
           logitRho = 0,
           beta0 = 3,
           gamma = rep(0,length(y)))

f<-function(par){
  getAll(data,par)
  rho = 2/(1 + exp(-logitRho))-1
  timeSteps <-length(gamma)
  sd_y <- exp(logsd[1])
  sd_gamma = exp(logsd[2])
  
  nll<-0
  
  nll <- nll -dnorm(gamma[1],0,sqrt(sd_gamma*sd_gamma/(1-rho*rho)),log=TRUE)
  for(i in 2:timeSteps){    
    nll <- nll -dnorm(gamma[i],rho*gamma[i-1],sd_gamma,log=TRUE)
  }

  logmu = beta0 + gamma
  
  nll = nll - sum(dnorm(y-exp(logmu), sd = sd_y,log = TRUE))
  
  ADREPORT(rho)
  ADREPORT(logmu)
  
  nll
}

