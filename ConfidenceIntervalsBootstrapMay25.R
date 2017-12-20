##########

#############COnfidence Intervals Bootstrappping

setwd("C:/Users/Palma/Desktop/MastersNewProject/Janbreak")

FinalCleanObese <- read.csv("FullObeseYear.csv")

#AGES

#YEARS

obseseprevalence<-function(datmat){
  y = as.integer(c(1988:2012))
  #########Automation of Code##############
  ##########Using Binomial Logit as link##########
  x<-matrix(c(23:90))   ##Create a matrix of ages(x)################
  #PropObese <-structure(c(FinalCleanObese$PropObese), .Dim = as.integer(c(68,25)))
  PropObese<-datmat
  xbar<-mean(x)  ##Create a matrix of deviations of ages(x) from the mean ages######### 
  ktguess<-matrix(0,3,25)
  FittedPropObeseAY<-matrix(0,68,25)
  Resi<-matrix(0,68,25)
  kt1guess<-0
  kt2guess<-0
  kt3guess<-0
  x_xbar2<-(x-xbar)^2
  #sigsq<-rep(var(x),length(x))
  #x_xbarsq<-((x-xbar)^2-sigsq)
  
  ##########Start Iteration to obtain coeffcients KT from the Generalised Linear Model Algorithm##########
  for (i in 1:25){  
    PropObeseAY <-PropObese[,i]
    #Obsesemod<-glm(PropObeseAY~x+x_xbar2,family=binomial(link ="logit"))
    Obsesemod<-lm(PropObeseAY~x+x_xbar2)   ##Check here again#
    kt1guess[i]<-Obsesemod$coefficients[1]
    kt2guess[i]<-Obsesemod$coefficients[2]
    kt3guess[i]<-Obsesemod$coefficients[3]
    FittedPropObeseAY[,i]<-cbind(Obsesemod$fitted.values)
  } 
  year<-1988:2012
  
  
  mx=xbar<-mean(x)  ##Create a matrix of deviations of ages(x) from the mean ages######### 
  mx2=mean((x-mx)^2)
  x_xbar2<-(x-xbar)^2-mx2
  
  
  
  n=length(x)  # number of ages
  m=length(y)  # number of years
  
  cy=(y[1]-x[n]):(y[m]-x[1])  # cohort approximate years of birth
  # initialise the parameter vectors
  wa=matrix(1,25,68)
  kap1v=(1:m)*0
  kap2v=kap1v
  kap3v=kap1v
  kap5v=kap1v    ###########New ###############
  beta2v=(x-mx)
  beta4v=(x-mx)^2-mx2
  beta3v=x*0+1
  gamma4v=(1:(n+m-1))*0
  
  
  ia=array((1:m),c(m,n))	 # matrix of year indexes, i, for the data
  ja=t(array((1:n),c(n,m)))  # matrix of age indexes, j, for the data
  ya=ia-ja		 	 # matrix of year of birth indexes for the data
  imj=(1-n):(m-1)		 # the range of values taken by i-j
  lg=n+m-1		 	 # number of different values taken by i-j
  ca=ya+y[1]-x[1]		 # matrix of years of birth
  
  
  # Now set weights to zero for cohorts with fewer than 5 observations
  for(k in 1:lg)
  {
    nk=sum((ca == cy[k])*wa)
    if(nk < 5)
    {
      wa=wa*(1- (ca == cy[k]))
    }
  }
  
  ww=cy*0+1  # this is a vector of 1's and 0's with
  # a 0 if the cohort is completely excluded
  for(k in 1:lg)
  {
    ww[k]=ww[k]*(sum((ca == cy[k])*wa) > 0)
  }
  ww1=order(-ww)[1]
  ww2=length(ww)+1-order(-ww[length(ww):1])[1]
  ww3=floor((ww1+ww2)/2)
  
  
  # Stage 1 
  # Full MLE, iterative
  iteration=0
  l0=-1000000
  l1=-999999
  # l1 is the latest estimate of the log-likelihood
  # l0 is the previous estimate
  # we continue to iterate if the improvement in log-likelihood
  # exceeds 0.0001
  fix.g3=0
  while(l1-l0 > 0.000000001)
  {
    iteration=iteration+1
    l0=l1
    
    if(fix.g3 == 0)
    {
      # Begin the calculations to apply the three constraints
      g3=gamma4v[ww1:ww2]
      cy3=cy[ww1:ww2]   ####Some prior  years 
      mm=mean(cy3+mx)
      s11=length(g3)
      s12=sum((cy3+mx-mm)^1)
      s13=sum((cy3+mx-mm)^2)
      s21=sum((cy3+mx-mm)^1)
      s22=sum((cy3+mx-mm)^2)
      s23=sum((cy3+mx-mm)^3)
      s31=sum((cy3+mx-mm)^2)
      s32=sum((cy3+mx-mm)^3)
      s33=sum((cy3+mx-mm)^4)
      t1=sum(g3)
      t2=sum(g3*(cy3+mx-mm))
      t3=sum(g3*(cy3+mx-mm)^2)
      ss=array(c(s11,s21,s31,s12,s22,s32,s13,s23,s33),c(3,3)) ### Creates a 3 by 3 matrix
      tt=array(c(t1,t2,t3),c(3,1))
      phi=(solve(ss) %*% tt)[,1]
      phi3=phi[1]-phi[2]*mm+phi[3]*mm^2
      phi4=phi[2]-2*phi[3]*mm
      phi5=phi[3]
      
      # apply the three constraints simultaneously
      gamma4v=gamma4v-phi3-phi4*(cy+mx)-phi5*((cy+mx)^2)
      kap1v=kap1v+phi3+phi4*y+phi5*(y^2-mx2)
      kap2v=kap2v-phi4-phi5*2*y
      kap3v=kap3v+phi5
    } # end the application of constraints if gamma4v is being estimated
    
    ##########Start Iteration to obtain coeffcients KT from the Linear Regression Model Algorithm##########
    for (i in 1:25){
      g3=gamma4v[(n+i-1):i]
      PropObeseAY <-PropObese[,i]-g3
      #Obsesemod<-glm(PropObeseAY~x+x_xbar2,family=binomial(link ="logit"))
      Obsesemod<-lm(PropObeseAY~beta2v+beta4v)
      kap1v[i]<-Obsesemod$coefficients[1]
      kap2v[i]<-Obsesemod$coefficients[2]
      kap3v[i]<-Obsesemod$coefficients[3]
      FittedPropObeseAY[,i]<-cbind(Obsesemod$fitted.values)
      Resi[,i]=cbind(Obsesemod$residuals)
    } 
    
    if(fix.g3 == 0)
    {
      for(k in 1:lg)
      {
        if(ww[k] == 0)
        { 
          gamma4v[k]=0 
        }	 # skip this k if the cohort is not included
        else 
        {
          id=0+(ya == imj[k])  # pick out entries with a specific year of birth
          
          gamma4v[k]=mean(Resi[id==1])
        }
      }
      
    }
    
    mhat=PropObese*0
    for(i in 1:m)
    {
      g3=gamma4v[(n+i-1):i]
      mhat[,i]=kap1v[i]+kap2v[i]*beta2v+kap3v[i]*beta4v+g3*beta3v
    }
    
    
    epsilon=PropObese-mhat  # matrix of  errors
    l1=sum(epsilon^2)   ## normal assumption estimate of  log likelihod
    cat(l1,"\n")
    
  }
  
  list(beta2v=beta2v,beta4v=beta4v,kap1v=kap1v,
       kap2v=kap2v, kap3v=kap3v,gamma4v=gamma4v,xv=x,yv=y,cy=cy,
       wa=wa,epsilon=epsilon,mhat=mhat,ll=l1)	
}

#install.packages("picante")
PropObese <-structure(c(FinalCleanObese$PropObese), .Dim = as.integer(c(68,25)))
fitobesemodel<-obseseprevalence(PropObese)
library(picante)
ktest<-matrix(NA,25,3)
nsim=1000  
KappMat=list()
BootResMat=list()   #######New
for(i in 1:nsim){ 
  randomiseError<-randomizeMatrix(fitobesemodel$epsilon,iterations = 1) # null.model=""richness"", frequency"
  bootfit <-fitobesemodel$mhat + randomiseError
  booftit1mod<-obseseprevalence(bootfit)
  BootResMat[[i]]<-list(booftit1mod$epsilon) #######New ##Should it be bootstap residuals Orignal vrs Fit or Boot Orig Vrs Fit
  ktest<-cbind(booftit1mod$kap1v,booftit1mod$kap2v,booftit1mod$kap3v)
  KappMat[[i]] <-list(ktest)
}

library(forecast)

#tsexp<-ts(BootResMat[[2]][[1]][20,],start=1988)
#plot(tsexp,type="b",col="blue")
#par(mfrow=c(1,1))
#acf(tsexp)
#pacf(tsexp)

#fit00<-Arima(tsexp,order=c(0,0,0))
#bestfor<-forecast(fit00,h=6)
#bestpred<-predict(fit00,h=6)   ######Output is just TS and SE  ## no TS predict 

#bestmod<-ar(tsexp,method ="mle")
#bestfor<-forecast(bestmod,h=6)
#plot(bestfor)



########Forecast the Bootstrap sample#############################

source("ForecastBootFunctions.R")
library(forecast)

kapvDataarray=list()
kap1vData.forecast<-rep(NA,6)
kap2vData.forecast<-rep(NA,6)
kap3vData.forecast<-rep(NA,6)
kapvSEaarray=list()
kap1vSE.fore<-rep(NA,6)
kap2vSE.fore<-rep(NA,6)
kap3vSE.fore<-rep(NA,6)

kapdatUCLarray=list()
kapdatLCLarray=list()
kap1vUCL<-rep(NA,6)
kap1vLCL<-rep(NA,6)
kap2vUCL<-rep(NA,6)
kap2vLCL<-rep(NA,6)
kap3vUCL<-rep(NA,6)
kap3vLCL<-rep(NA,6)

resDataarray=list()

for(i in 1:nsim)
{
  
  ##########################Obtain Kap Predictions
  kap1vData.forecast<-sarima.forboot(KappMat[[i]][[1]][,1],6,0,1,0)$pred[1:6]
  kap2vData.forecast<-sarima.forboot(KappMat[[i]][[1]][,2],6,0,1,0)$pred[1:6]
  kap3vData.forecast<-sarima.forboot(KappMat[[i]][[1]][,3],6,0,1,0)$pred[1:6]
  kapfor<-rbind(kap1vData.forecast,kap2vData.forecast, kap3vData.forecast)
  kapvDataarray[[i]] <-list( kapfor)
  
  ##########################Obtain Standard Errors
  kap1vSE.fore<-sarima.forboot(KappMat[[i]][[1]][,1],6,0,1,0)$se[1:6]
  kap2vSE.fore<-sarima.forboot(KappMat[[i]][[1]][,2],6,0,1,0)$se[1:6]
  kap3vSE.fore<-sarima.forboot(KappMat[[i]][[1]][,3],6,0,1,0)$se[1:6]
  kapforInt<-rbind(kap1vSE.fore,kap2vSE.fore, kap3vSE.fore)
  kapvSEaarray[[i]] <-list(kapforInt)
  
  #######################Obtain Upper Limits 
  kap1vUCL<-kap1vData.forecast+ 2*kap1vSE.fore 
  kap2vUCL<-kap2vData.forecast+ 2*kap2vSE.fore 
  kap3vUCL<-kap3vData.forecast+ 2*kap3vSE.fore 
  kapdatUCL<-rbind(kap1vUCL,kap2vUCL,kap3vUCL)
  kapdatUCLarray[[i]]<-list(kapdatUCL)
  
  #######################Obtain Lower Limits
  
  kap1vLCL<-kap1vData.forecast - 2*kap1vSE.fore 
  kap2vLCL<-kap2vData.forecast- 2*kap2vSE.fore 
  kap3vLCL<-kap3vData.forecast- 2*kap3vSE.fore 
  kapdatLCL<-rbind(kap1vLCL,kap2vLCL,kap3vLCL)
  kapdatLCLarray[[i]]<-list( kapdatLCL)
  
  #######################TS Model for residuals  Lower Limits
  restsmod <-sarima.forboot(BootResMat[[i]][[1]][1,],6,0,0,0)$pred[1:6]
  resse <- sarima.forboot(BootResMat[[i]][[1]][1,],6,0,0,0)$se[1:6]
  resseUL<-restsmod+2*resse
  resseLL<-restsmod-2*resse
  respred<-rbind(restsmod,resseUL,resseLL)
  
  ##########################Try using forecast package ##########
  # restsmodts<-ts(BootResMat[[i]][[1]][1,], start=1988,frequency=1)
  #restsmod <-ar(restsmodts,method="mle")
  ###restsmodfor<-predict(restsmod,h=6)  Output is Prediction and SE not Prediction intervals
  # resfor<- forecast(restsmod ,h=6)
  #respred<-rbind(resfor$mean[1:6],resfor$upper[,2],resfor$lower[,2]) ##95% interval
  resDataarray[[i]] <-list(respred)
}  


###############################  Fit the forecast ####################
mhatforarray=list()
mhatfor<-matrix(0,68,6)
for(i in 1:nsim){ 
  for(j in 1:6){
    mhatfor[,j] <- kapvDataarray[[i]][[1]][1,j]+kapvDataarray[[i]][[1]][2,j]*fitobesemodel$beta2v+kapvDataarray[[i]][[1]][3,j]* fitobesemodel$beta4v
    +  resDataarray[[i]][[1]][1,j]
    mhatforarray[[i]] <-list(mhatfor)
  }
}


###############################  Fit the Uper Pred Limits ####################
mhatUCLarray=list()
mhatUCLfor<-matrix(0,68,6)
for(i in 1:nsim){ 
  for(j in 1:6){
    mhatUCLfor[,j] <- kapdatUCLarray[[i]][[1]][1,j]+kapdatUCLarray[[i]][[1]][2,j]*fitobesemodel$beta2v+kapdatUCLarray[[i]][[1]][3,j]* fitobesemodel$beta4v
    + resDataarray[[i]][[1]][2,j]
    mhatUCLarray[[i]] <-list( mhatUCLfor)
  }
}

###############################  Fit the lOWER Pred Limits ####################
mhatLCLarray=list()
mhatLCLfor<-matrix(0,68,6)
for(i in 1:nsim){ 
  for(j in 1:6){
    mhatLCLfor[,j] <- kapdatLCLarray[[i]][[1]][1,j]+kapdatLCLarray[[i]][[1]][2,j]*fitobesemodel$beta2v+kapdatLCLarray[[i]][[1]][3,j]* fitobesemodel$beta4v
    + resDataarray[[i]][[1]][3,j]
    mhatLCLarray[[i]] <-list( mhatLCLfor)
  }
}


################ Plot Forecasts 
###Figure 6 

par(mfrow=c(2,2))


plot(1988:2012,fitobesemodel$mhat[23,] ,col="black", pch=20, type="l",ylim=c(.10,.39), xlim=c(1988,2020),main
     ="6 Year forecast for Age 45",xlab="Time", ylab="Proportion Obese")
points( 1988:2012,fitobesemodel$mhat[23, ],col="black",pch=1)

lines( 2013:2018,mhatforarray[[500]][[1]][23,],col="red",lwd=2)
points(2013:2018,mhatforarray[[500]][[1]][23,],col="red",pch=20)

lines( 2013:2018,mhatLCLarray[[25]][[1]][23,],col="blue",lwd=2)
points(2013:2018,mhatLCLarray[[25]][[1]][23,],col="blue",pch=1)

lines( 2013:2018,mhatUCLarray[[975]][[1]][23,],col="blue",lwd=2)
points(2013:2018,mhatUCLarray[[975]][[1]][23,],col="blue",pch=1)



#############################################
#Plot for age 55
plot(1988:2012,fitobesemodel$mhat[33,] ,col="black", pch=20, type="l",ylim=c(.10,.43), xlim=c(1988,2020),main
     ="6 Year forecast for Age 55",xlab="Time", ylab="Proportion Obese")
points( 1988:2012,fitobesemodel$mhat[33,],col="black",pch=1)


lines( 2013:2018,mhatforarray[[500]][[1]][33,],col="red",lwd=2)
points(2013:2018,mhatforarray[[500]][[1]][33,],col="red",pch=20)

lines( 2013:2018,mhatLCLarray[[25]][[1]][33,],col="blue",lwd=2)
points(2013:2018,mhatLCLarray[[25]][[1]][33,],col="blue",pch=1)

lines( 2013:2018,mhatUCLarray[[975]][[1]][33,],col="blue",lwd=2)
points(2013:2018,mhatUCLarray[[975]][[1]][33,],col="blue",pch=1)



#Plot for age 65
plot(1988:2012,fitobesemodel$mhat[43,] ,col="black", pch=20, type="l",ylim=c(.10,.42), xlim=c(1988,2020),main
     ="6 Year forecast for Age 65",xlab="Time", ylab="Proportion Obese")
points(1988:2012,fitobesemodel$mhat[43,],col="black",pch=1)

lines( 2013:2018,mhatforarray[[500]][[1]][43,],col="red",lwd=2)
points(2013:2018,mhatforarray[[500]][[1]][43,],col="red",pch=20)

lines( 2013:2018,mhatLCLarray[[25]][[1]][43,],col="blue",lwd=2)
points(2013:2018,mhatLCLarray[[25]][[1]][43,],col="blue",pch=1)

lines( 2013:2018,mhatUCLarray[[975]][[1]][43,],col="blue",lwd=2)
points(2013:2018,mhatUCLarray[[975]][[1]][43,],col="blue",pch=1)





#Plot for age 75
plot(1988:2012,fitobesemodel$mhat[53,] ,col="black", pch=20, type="l",ylim=c(.10,.38), xlim=c(1988,2020),main
     ="6 Year forecast for Age 75",xlab="Time", ylab="Proportion Obese")
points(1988:2012,fitobesemodel$mhat[53,],col="black",pch=1)

lines( 2013:2018,mhatforarray[[500]][[1]][53,],col="red",lwd=2)
points(2013:2018,mhatforarray[[500]][[1]][53,],col="red",pch=20)

lines( 2013:2018,mhatLCLarray[[25]][[1]][53,],col="blue",lwd=2)
points(2013:2018,mhatLCLarray[[25]][[1]][53,],col="blue",pch=1)

lines( 2013:2018,mhatUCLarray[[975]][[1]][53,],col="blue",lwd=2)
points(2013:2018,mhatUCLarray[[975]][[1]][53,],col="blue",pch=1)




