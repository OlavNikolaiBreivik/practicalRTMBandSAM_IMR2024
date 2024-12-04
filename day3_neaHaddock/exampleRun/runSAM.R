library(stockassessment)

cn<-read.ices("exampleRun/cn.dat") #Catch numbers
cw<-read.ices("exampleRun/cw.dat") #Catch weight
dw<-read.ices("exampleRun/dw.dat") #Discard weight
lf<-read.ices("exampleRun/lf.dat") #Landing fraction
lw<-read.ices("exampleRun/lw.dat") #Landing weight
mo<-read.ices("exampleRun/mo.dat") #Maturity
nm<-read.ices("exampleRun/nm.dat") #Natural mortality
pf<-read.ices("exampleRun/pf.dat") #Proportion fishing mortality before calculating SSB
pm<-read.ices("exampleRun/pm.dat") #Proportion natural mortality before calculating SSB
sw<-read.ices("exampleRun/sw.dat") #stock weight
surveys<-read.ices("exampleRun/survey.dat") #Surveys

dat<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn, 
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)


conf = loadConf(dat,file= "exampleRun/conf.cfg")
par = defpar(dat,conf)
fit = sam.fit(dat,conf,par)


##### include outlier in given year, fleet and age####
year = 1970
fleet = 1
age = 7

aux = as.data.frame(fit$data$aux) #data frame with years etc. ordered in same way as observations

dataOutlier = fit$data
obsToModify = which(aux$year==year& aux$fleet==fleet& aux$age == age)
dataOutlier$logobs[obsToModify] = dataOutlier$logobs[obsToModify] +10  #Include an outlier in given year, age and fleet
conf = loadConf(dat,file= "exampleRun/conf.cfg")
par = defpar(dataOutlier,conf)
fitOutlier = sam.fit(dataOutlier,conf, par)

ssbplot(c(fitOfficial = fit,fitOutlier = fitOutlier),addCI = TRUE)
  
confRobust = conf
confRobust$fracMixObs[1] = 0.01
fitOutlierRobust = sam.fit(dataOutlier,confRobust, par)

ssbplot(c(fitOfficial = fit,fitOutlier = fitOutlier,fitRobust = fitOutlierRobust),addCI = TRUE)

