N0 <- 1000; F <- 0.4; M <- 0.2
C0 <- F/(F+M)*(1-exp(-F-M))*N0
N1 <- N0*exp(-F-M)

C1 <- F/(F+M)*(1-exp(-F-M))*N1
N2 <- N1*exp(-F-M)

N1
N2




#Solution:

#Do it in a for-loop
N0<-1000; F<-.4; M<-.2
nYears = 10
N<-c(N0, rep(NA,nYears))
C = rep(0,nYears)

for(i in 1:(length(N)-1)){
  C[i] <- F/(F+M)*(1-exp(-F-M))*N[i]
  N[i+1] = N[i]*exp(-F-M)
}

round(N)
round(C)


#Percentage that dies from fishery in this scenario
(1-exp(-F-M)) - (1-exp(-M))
