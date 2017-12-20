########################################################################
####  Iteration to fit the Obsesity data 
########################################################################


FinalCleanObese <- read.csv("FullObeseYear.csv")

#AGES

#YEARS
y = as.integer(c(1988:2012))

#########Automation of Code##############
##########Using Binomial Logit as link##########
x<-matrix(c(23:90))   ##Create a matrix of ages(x)################
PropObese <-structure(c(FinalCleanObese$PropObese), .Dim = as.integer(c(68,25)))
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
  if(nk < 1)
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
  
  
  epsilon=PropObese-mhat  # matrix of standardised errors
  l1=sum(epsilon^2)   ## normal assumption estimate of  log likelihod
  cat(l1,"\n")
  
}


year<-1988:2012
ktguess<-cbind(year,kap1v,kap2v,kap3v)
#ktguess
#FittedPropObeseAY
#View(epsilon)


####Graph of results

#########Figure 6   
par(mfrow=c(1,1))
with(FinalCleanObese,plot(ObeseCount/Number~AGE, pch =20,col="gray",
                          main="Proposed Model Fit ",
                          ylim=c(0,.35),bty="l",xlim=c(23,94),
                          xlab="AGE",
                          ylab="Proportion Obese"))

for (i in 1:m){ 
  lines(23:90,mhat[,i],type="l",
        xlim=c(23,94),ylim=c(0,.35),bty="l",col="black",lwd= 1.7)
}
grid(nx=NULL,ny=NULL)
text(x=93.5, y=seq(0.02411201,0.1105592,.02), labels=c(1988,1993,2000,2005,2012)) 



#max(MedPloTestFit[68,])
#min(MedPloTestFit[68,])
#plot(x,PropObese[,15],pch=20)
#lines(x,mhat[,15],col="blue",lwd=3)

MAPENEWMODEL<-mean(abs(mhat-PropObese)/PropObese)
MAPENEWMODEL


MAPENEWMODELana<-mean(abs(AnovaFittedPropObeseAY-PropObese)/PropObese)
MAPENEWMODELana

npar=length(kap1v)+length(kap2v)+length(kap3v)+(1-fix.g3)*(length(gamma4v)-3)



###Proposed model##########
sig2<-mean(epsilon^2)

npar=length(kap1v)+length(kap2v)+length(kap3v)+(1-fix.g3)*(length(gamma4v)-3)
BIC = -2*log(sig2)+ npar*log(n*m)
BIC
AIC =-2*log(sig2)+2*npar
AIC






############residuals Vrs Fitted



###########Plot of coefficients###################
###Figure 3 
par(mfrow=c(1,3))
plot(y,kap1v,xlab="Year",ylab="kappa1", type="o",main="Parameter Kt(1)",lwd=1,col="black",pch=20)
grid(nx=NULL,ny=NULL)
plot(y,kap2v,xlab="Year",ylab="kappa2", type="o",main="Parameter Kt(2)",lwd=1,col="black",pch=20)
grid(nx=NULL,ny=NULL)
plot(y,kap3v,xlab="Year",ylab="kappa3", type="o",main="Parameter Kt(3)",lwd=1,col="black",pch=20)
#FinalCleanObese$cohyr<-FinalCleanObese$IYEAR-FinalCleanObese$AGE
grid(nx=NULL,ny=NULL)

############Cohort Effect Plot#########################
##Figure 4
par(mfrow=c(1,1))
plot(1898:1989,gamma4v,xlab="Cohort Years",ylab="gamma",type="o",pch=20,lwd=1,col="black",main="Cohort Effect Estimates")
grid(nx=NULL,ny=NULL)



######################Normal Regression without cohort Fit#########
FittedPropObeseAYREG<-matrix(0,68,25)
kt1guess<-0
kt2guess<-0
kt3guess<-0
beta2v=(x-mx)
beta4v=(x-mx)^2-mx2
for (i in 1:25){  
  PropObeseAY <-PropObese[,i]
  ObsesemodREG<-lm(PropObeseAY~beta2v+beta4v)   ##Check here again#
  kt1guess[i]<-ObsesemodREG$coefficients[1]
  kt2guess[i]<-ObsesemodREG$coefficients[2]
  kt3guess[i]<-ObsesemodREG$coefficients[3]
  FittedPropObeseAYREG[,i]<-cbind(ObsesemodREG$fitted.values)
} 

#############Figure 5
par(mfrow=c(1,1))
with(FinalCleanObese,plot(ObeseCount/Number~AGE,pch=20,col="gray",
                          main="Quadratic Regression with no Cohort Fit ",
                          ylim=c(0,.35),bty="l", xlim=c(23,94),
                          xlab="AGE",
                          ylab="Proportion Obese"))

for (i in 1:25){ 
  lines(23:90,FittedPropObeseAYREG[,i],type="l",
        xlim=c(23,94), ylim=c(0,.35),bty="l",col="black",lwd=1.7)
}
grid(nx=NULL,ny=NULL)
text(x=93.5, y=seq(0.02465707, 0.1087422,.02), labels=c(1988,1993,2000,2005,2012)) 

max(FittedPropObeseAYREG[68,])

min(FittedPropObeseAYREG[68,])
########################################################################


##########Figure 1 
par(mfrow=c(1,1))
with(FinalCleanObese,plot(ObeseCount/Number~AGE,pch=20,col="gray",
                          main="Obese Proportions for the U.S. by Age from 1988 to 2012",
                          ylim=c(0,.35),bty="l", xlim=c(23,94),
                          xlab="AGE",
                          ylab="Proportion Obese"))

for (i in 1:25){ 
  lines(23:90,FittedPropObeseAYREG[,i],type="l",
        xlim=c(23,94), ylim=c(0,.35),bty="l",col="black",lwd=1.7)
}
grid(nx=NULL,ny=NULL)
text(x=93.5, y=seq(0.02465707, 0.1087422,.02), labels=c(1988,1993,2000,2005,2012)) 

###################################################

####Figure 2

with(FinalCleanObese,plot(ObeseCount/Number~IYEAR,pch=20,col="gray",
                          main="Obese Proportions for the U.S. from 1988 to 2012 by Age",xlim=c(1988,2014),
                          ylim=c(0,.35),bty="l",
                          xlab="YEAR",
                          ylab="Proportion Obese"))

for (i in 1:68){ 
  lines(1988:2012,FittedPropObeseAYREG[i,],type="l",
        xlim=c(1988,2014), ylim=c(0,.35),bty="l",col="black",lwd=1.5)
}
grid(nx=NULL,ny=NULL)

text(x=2013.5, y=seq(0.1105592,0.3397288,.05), labels=c(90,85,50,40,35)) 

##########################Graph by Year ######################

######Figure 7

with(FinalCleanObese,plot(ObeseCount/Number~IYEAR,pch=20,col="gray",
                          main="Proposed Model Fit",xlim=c(1988,2014),
                          ylim=c(0,.35),bty="l",
                          xlab="YEAR",
                          ylab="Proportion Obese"))

for (i in 1:68){ 
  lines(1988:2012,mhat[i,],type="l",
        xlim=c(1988,2014), ylim=c(0,.35),bty="l",col="black",lwd=1.5)
}
grid(nx=NULL,ny=NULL)

text(x=2013.5, y=seq( 0.1105592,0.3397288,.05), labels=c(90,85,50,40,35)) 

max(mhat[,25])                            

min(mhat[,25])



### ANOVAConstraint Based Approach ##########

FinalCleanObese <- read.csv("FullObeseYear.csv")
###########################################################################
#########Automation of Code########################
FinalCleanObese$cohort<-FinalCleanObese$IYEAR-FinalCleanObese$AGE

Obsesemod.lm<-lm(log(PropObese)~factor(AGE)+factor(IYEAR)+factor(cohort),data=FinalCleanObese)

summary(Obsesemod.lm)


AnovaFittedPropObeseAY <-structure(c(exp(Obsesemod.lm$fitted.values)), .Dim = as.integer(c(68,25)))

################Figure 8

with(FinalCleanObese,plot(ObeseCount/Number~AGE,pch=20,col="gray",
                          main="Constraint Based AGE-PERIOD-COHORT Fit ",
                          ylim=c(0,.35),bty="l",xlim=c(23,94),
                          xlab="AGE",
                          ylab="Proportion Obese"))

for (i in 1:25){ 
  lines(23:90,AnovaFittedPropObeseAY[,i],type="l",
        xlim=c(23,94), ylim=c(0,.35),bty="l",col="black",lwd=1.7)
}
grid(nx=NULL,ny=NULL)
text(x=93.5, y=seq( 0.01876361,0.09034112,.02), labels=c(1988,1993,2000,2012)) 

########################################################


#######################Median Polish 
#install.packages("STMedianPolish")
library(STMedianPolish)

FinalCleanObese <- read.csv("FullObeseYear.csv")

PropObese <-structure(c(FinalCleanObese$PropObese), .Dim = as.integer(c(68,25)))
logPropObese<-log10(PropObese)
MedPloTest<- MedianPolishM(logPropObese)

MedPloTestFit<-(logPropObese-MedPloTest$residuals)


MAPENEWMODEMedpol<-mean(abs(MedPloTestFit-logPropObese))
MAPENEWMODEMedpol


l1medpol<-mean(MedPloTest$residuals^2)
npar=length(MedPloTest$effects[[1]])+length(MedPloTest$effects[[2]])+ 1
npar
BIC = -2*log(l1medpol)+ npar*log(n*m)
BIC
AIC =-2*log(l1medpol)+2*npar
AIC


###Proposed model##########
sig2<-mean(epsilon^2)

npar=length(kap1v)+length(kap2v)+length(kap3v)+(1-fix.g3)*(length(gamma4v)-3)
BIC = -2*log(sig2)+ npar*log(n*m)
BIC
AIC =-2*log(sig2)+2*npar
AIC


########Plot Data and fitted#############

#############Figure 9
par(mfrow=c(1,1))
with(FinalCleanObese,plot(ObeseCount/Number~AGE,pch=20,col="gray",
                          main="MEDIAN POLISH Model Fit ",
                          ylim=c(0,.35),bty="l",xlim=c(23,94),
                          xlab="AGE",
                          ylab="Proportion Obese"))

for (i in 1:25){ 
  lines(23:90,MedPloTestFit[,i],type="l",
        xlim=c(23,94),ylim=c(0,.35),bty="l",col="black",lwd=1.7)
}
grid(nx=NULL,ny=NULL)
text(x=93.5, y=seq( -0.03164201,0.1414252,.04), labels=c(1988,1993,2000,2005,2012)) 



##########################Time Series Analysis##################
# Select  the Age group 45 #############

###########################################################################################
#Time Series analysis of coefficients #################

library(forecast)
source("ForecastDrKimFunctions.R")  ###Dr Kim forecast Package ##########
source("ForecastNewFunctions.R")


kt1<-kap1v
kt2<-kap2v
kt3<-kap3v
gamma4<-gamma4v

################################################################

#########Figure 10
par(mfrow=c(1,3))
sarima.forNew(kt1,6,0,1,0,datstart=1988:2018,xstend=2013:2018)
sarima.forNew(kt2,6,0,1,0,datstart=1988:2018,xstend=2013:2018)
sarima.forNew(kt3,6,0,1,0,datstart=1988:2018,xstend=2013:2018)
sarima.forNew(gamma4,6,0,0,0,datstart=1898:1995,xstend=1991:1996)

arima(gamma4)

gamma4<-ts(gamma4)


###################################################################33


###############Figure 11  ### Found in COnfidence Interval boostrap R file
###################################################################


