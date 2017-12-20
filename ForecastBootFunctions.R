

sarimaboot=function(data,p,d,q,P=0,D=0,Q=0,S=-1){ 
  n=length(data)
  constant=1:n
  xmean=matrix(1,n,1)
  if (d>0)  
    fitit=arima(data, order=c(p,d,q), seasonal=list(order=c(P,D,Q), period=S),xreg=constant,include.mean=F)
  if (d<.00001)
    fitit=arima(data, order=c(p,d,q), seasonal=list(order=c(P,D,Q), period=S),xreg=xmean,include.mean=F)
  if (d+D>1)
    fitit=arima(data, order=c(p,d,q), seasonal=list(order=c(P,D,Q), period=S))
  if (S < 0) goof=20 else goof=3*S
  tsdiag(fitit,gof.lag=goof)
  k=length(fitit$coef)
  BIC=log(fitit$sigma2)+(k*log(n)/n)
  AICc=log(fitit$sigma2)+((n+k)/(n-k-2))
  AIC=log(fitit$sigma2)+((n+2*k)/n)
  list(fit=fitit, AIC=AIC, AICc=AICc, BIC=BIC)
}

sarima.forboot=function(xdata,nahead,p,d,q,P=0,D=0,Q=0,S=-1,tol=.001){ 
  data=as.ts(xdata) 
  n=length(data)
  constant=  1:n  ##1:n
  xmean=matrix(1,n,1)
  if (d>0 & D>0) 
    fitit=arima(data, order=c(p,d,q), seasonal=list(order=c(P,D,Q), period=S))
  if (d>0 & D==0)  
    fitit=arima(data, order=c(p,d,q), seasonal=list(order=c(P,D,Q), period=S),
                xreg=constant,include.mean=F)
  if (d==0 & D==0)
    fitit=arima(data, order=c(p,d,q), seasonal=list(order=c(P,D,Q), period=S),
                xreg=xmean,include.mean=F)
  if (d==0 & D>0)  
    fitit=arima(data, order=c(p,d,q), seasonal=list(order=c(P,D,Q), period=S),
                xreg=constant,include.mean=F)
  if (d>0 & D>0)   nureg=NULL
  if (d>0 & D==0)  nureg=(n+1):(n+nahead)
  if (d==0 & D==0) nureg=matrix(1,nahead,1)
  if (d==0 & D>0)  nureg=(n+1):(n+nahead)
  fore=predict(fitit, n.ahead=nahead, newxreg=nureg)  
  #-- graph:
  U = fore$pred + 2*fore$se
  L = fore$pred - 2*fore$se
  a=max(1,n-100)
  minx=min(data[a:n],L)
  maxx=max(data[a:n],U)
  #t1=xy.coords(data)$x; 
  t1=1:length(data);
  if(length(t1)<101) strt=t1[1] else strt=t1[length(t1)-100]
  #t2=xy.coords(fore$pred)$x; 
  t2=length(data)+c(1:nahead)
  endd=t2[length(t2)]
  xllim=c(strt,endd)
 # ts.plot(data,fore$pred,col=1:2, xlim=xllim, ylim=c(minx,maxx), ylab=deparse(substitute(xdata))) 
  #lines(fore$pred, col="red", type="p")
 # lines(U, col="blue", lty="dashed")
  #lines(L, col="blue", lty="dashed")
  #
  return(fore)
  return( U)
  return(L)
}