fitWeb = stockassessment::fitfromweb("NEACod_2024_Final")
fitWeb = stockassessment:::refit(fitWeb)


# function for cross-validation  
xval <- function(fit, year=NULL, fleet=NULL, age=NULL, ...){
  data <- fit$data
  nam <- c("year", "fleet", "age")[c(length(year)>0,length(fleet)>0,length(age)>0)]
  if((length(year)==0) & (length(fleet)==0) & (length(age)==0)){
    idx <- rep(TRUE,nrow(data$aux))
  }else{
    idx <- !do.call(paste, as.data.frame(data$aux[,nam,drop=FALSE])) %in% do.call(paste, as.data.frame(cbind(year=year, fleet=fleet, age=age)))
  }
  idx <- !idx
  data$logobs[idx] <- NA
  idx2 <- which(is.na(data$logobs))
  conf <- fit$conf
  par <- defpar(data,conf)
  thisfit <- sam.fit(data, conf, par, rm.unidentified = TRUE, silent=TRUE,...)
  ret <- as.data.frame(cbind(data$aux[idx2,], obs=fit$data$logobs[idx2], pred=thisfit$pl$missing, predSd=thisfit$plsd$missing))
  ret <- ret[complete.cases(ret),]
  attr(ret, "fit") <- thisfit
  return(ret)
}



pred = xval(fitWeb,year = c(2000:2005,2023:2024))

plot(pred$obs,pred$pred)
abline(0,1)

png("crossValSSB.png")
ssbplot(fitWeb,xlim = c(1990,2024), main = "SSB COD",cicol = rgb(0,0, 0, alpha = 0.25),cex.main = 2, cex.lab = 1.4)
ssbplot(attributes(pred)$fit, add = TRUE, cicol =rgb(0, 0, 1, alpha = 0.25),col ='blue')
dev.off()

png("crossVal.png")
plot(pred$obs,pred$pred, main = "Prediction vs Truth")
abline(0,1)
dev.off()
