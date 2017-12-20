
sarima.forNew=function(xdata,nahead,p,d,q,P=0,D=0,Q=0,S=-1,tol=.001,datstart,xstend){ 
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
  
  ndata<-c(xdata,fore$pred[1:nahead])
  
  plot(datstart,ndata, ylim=c(minx,maxx), xlab="Time", ylab=deparse(substitute(xdata)),pch=20,type="l") 
  lines(xstend,fore$pred[1:nahead], col="red", type="p",pch=20)
  lines(xstend,c(U[1:nahead]), col="blue", lty="dashed")
  lines(xstend,c(L[1:nahead]), col="blue", lty="dashed")
  grid(nx=NULL,ny=NULL)

  return(fore)  

}


