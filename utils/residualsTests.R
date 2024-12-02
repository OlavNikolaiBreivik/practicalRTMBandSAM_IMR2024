##' Residual diagnostic tests for a SAM model fit. Function created by Casper W. Berg, thank you for sharing. 
##'
##' 
##' @title Residual diagnostic tests for a SAM model fit
##' @param fit an object of class ‘sam’ 
##' @param resid residuals (object of class "samres") from 'fit'. Calculated from 'fit' if not provided here.  
##' @param bubbles show bubble plots (default: TRUE)
##' @param palette vector of colors used for p-value groupings
##' @param p.breaks breaks for p-value groupings
##' @param reverseFleetOrder reverse order of fleets within year for OSA residuals (default: FALSE) 
##' @return list of pvalues
##' @details ...
residDiagPlot<-function(fit,resid=NULL,bubbles=TRUE,palette=hcl.colors(4,palette="Blues"),p.breaks=c(0,0.001,0.05,0.1,1),reverseFleetOrder=FALSE){
  
  stopifnot( (length(palette)+1)==length(p.breaks))
  if(is.null(resid)){
    if(reverseFleetOrder){
      resid = residuals(fit,subset=order(fit$data$aux[,"year"],max(fit$data$aux[,"fleet"])-fit$data$aux[,"fleet"],fit$data$aux[,"age"]))
    } else {
      resid = residuals(fit)
    }
  }
  resid.df = data.frame(residual=resid$residual,fleet=resid$fleet,age=resid$age,year=resid$year)
  resid.df$residual[ is.nan(resid.df$residual) ] = NA
  
  ## Biomass indices have age = -1, change to min age for convenience
  if(any(resid.df$age<0)){
    resid.df$age[ resid.df$age<0 ] = fit$conf$minAge 
  }    
  restab = xtabs(residual ~ year + age + fleet,data=resid.df)
  restab[ restab==0 ] = NA
  
  fnames <- sapply(attr(fit$data,"fleetNames"),substr,start=0,stop=14)
  
  acf.time.p <- bias <- variance <- bias.p <- variance.p <- matrix(NA,nrow=dim(restab)[3],ncol=fit$conf$maxAge - fit$conf$minAge + 1,dimnames=list(Fleet=fnames,Age=fit$conf$minAge:fit$conf$maxAge))
  
  shapiro <- acf.age.p <- meanvar <- meanvar.p <- numeric(dim(restab)[3])
  
  chisqtest<-function(x,testvar=1){
    n <- length(x)
    s2 <- var(x)
    chsq <- (n-1)*s2/testvar
    pval <- 1 - pchisq(chsq,n-1)
    pval
  }
  
  noagesf <- fit$data$maxAgePerFleet - fit$data$minAgePerFleet
  
  bp<-function(x,y,v, scale=3, ...){
    xlim <- c(min(x) - 1, max(x) + 1)
    ylim <- c(min(y) - 1, max(y) + 1)
    plot(x,y,cex=sqrt(abs(v))*scale, col=ifelse(v<0,rgb(1, 0, 0, alpha = 0.5),rgb(0, 0, 1, alpha = 0.5)), pch=19, xlim=xlim, ylim=ylim,...)
    points(x[v>0],y[v>0],cex=sqrt(v[v>0])*scale, col=rgb(0, 0, 1, alpha = 0.5), pch=19, ...)
  }
  
  par(mfrow=n2mfrow(dim(restab)[3]+6),las=1,mar=c(5,8,5,5))
  
  for(fl in 1:dim(restab)[3]){
    ## fleet based tests
    
    flsel = resid$fleet==fl
    if(bubbles){
      bp(resid$year[flsel],resid$age[flsel],resid$residual[flsel],main=paste(attr(fit$data,"fleetNames")[fl],"(",fl,")"),xlab="Year",ylab="Age")
    }
    shapiro[fl] <- shapiro.test(na.omit(as.vector(restab[,,fl])))$p.value        
    
    acfvec <- as.vector(apply(restab[,,fl],1,function(x) as.numeric(c(as.vector(x),rep(NA,noagesf[fl]+1)) )))
    
    
    if(is.list(acfvec)) acfvec <- do.call("c",acfvec)
    acf.age.p[fl] <- Box.test( acfvec,type="Ljung-Box",lag=1)$p.value
    
    sel <- resid.df$fleet==fl
    sel2 <- fit$rep$predObs[sel] < median(fit$rep$predObs[sel])
    lofit <- resid.df$residual[sel][ sel2 ]
    hifit <- resid.df$residual[sel][ !sel2 ]
    
    meanvar[fl] <- var(hifit,na.rm=TRUE) / var(lofit,na.rm=TRUE)
    meanvar.p[fl] <- pf( meanvar[fl], length(hifit)-1, length(lofit)-1)
    
    for(a in fit$data$minAgePerFleet[fl]:fit$data$maxAgePerFleet[fl]){
      ## fleet and age based tests
      colu = ifelse(fit$data$fleetTypes[fl]==3,1,a-fit$conf$minAge+1) ## take care of biomass indices (no ages)
      bias[fl,colu] = mean(restab[,colu,fl],na.rm=TRUE)
      variance[fl,colu] = var(na.omit(restab[,colu,fl]))
      bias.p[fl,colu] = t.test(restab[,colu,fl])$p.value
      variance.p[fl,colu] = chisqtest(na.omit(restab[,colu,fl]))
      acf.time.p[fl,colu] = Box.test( restab[,colu,fl],type="Ljung-Box")$p.value
      ## Note, that this is also testing for negative autocorrelation.
      ## Consider only testing for positive? (or use e.g. red colors instead of blue for significant negative correlations?) 
    }
    
  }
  
  plot.matrix:::plot.matrix(apply(t(bias.p),2,rev),breaks=p.breaks,col=palette,fmt.cell = "%.4f",xlab="Fleet",ylab="Age",main="Bias")
  plot.matrix:::plot.matrix(apply(t(variance.p),2,rev),breaks=p.breaks,col=palette,fmt.cell = "%.4f",xlab="Fleet",ylab="Age",main="Variance")
  
  plot.matrix:::plot.matrix(apply(t(acf.time.p),2,rev),breaks=p.breaks,col=palette,fmt.cell = "%.4f",xlab="Fleet",ylab="Age",main="ACF time direction")
  
  plot.matrix:::plot.matrix(t(matrix(shapiro,nrow=1)),breaks=p.breaks,col=palette,fmt.cell = "%.4f",xlab="",ylab="fleet",main="Normality",key=NULL)
  
  plot.matrix:::plot.matrix(t(matrix(acf.age.p,nrow=1)),breaks=p.breaks,col=palette,fmt.cell = "%.4f",xlab="",ylab="fleet",main="Correlation age direction",key=NULL)
  
  plot.matrix:::plot.matrix(t(matrix(meanvar.p,nrow=1)),breaks=p.breaks,col=palette,fmt.cell = "%.4f",xlab="",ylab="fleet",main="Mean-variance relationship",key=NULL)
  
  list(bias.p,variance.p,acf.time.p,shapiro,acf.age.p,meanvar.p)
  
}


