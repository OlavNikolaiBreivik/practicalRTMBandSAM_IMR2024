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


plot(fit)
ssbplot(fit)
fbarplot(fit)
recplot(fit)
fselectivityplot(fit)


