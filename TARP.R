# Here we compute RIS-RP and RIS-PCR for normal-linear model

######################### use ##############################
#@input: x : standardized taining design matrix in nXp format
#@input: y : training response, a n-vector
#@input: x_test : standardized test design matrix in nXp format

RISRP_res=RIS_PCR(x,y,x_test,alpha=0.95)
y_hat=RIS_PCR(x,y,x_test,alpha=0.95)
y_hat=RISRP_res[,1]     # predicted value
lower_PI=RISRP_res[,2]  # lower 100(1-alpha)% confidence limit
upper_PI=RISRP_res[,3]  # upper 100(1-alpha)% confidence limit


plot(y_test,type="l",lwd=3,ylim=c(-10,10),xlab="Predicted and observed values")
par(new=TRUE)
plot(y_hat,type="l",lwd=3,ylim=c(-10,10),col="blue")
par(new=TRUE)
plot(lower_PI,type="l",lwd=3,ylim=c(-10,10),col="red",lty=2)
par(new=TRUE)
plot(upper_PI,type="l",lwd=3,ylim=c(-10,10),col="red",lty=2)

################### Functions ###################
rRM.MBCR=function(k,psi)
{
	p1=psi;	p2=(1-2*psi); p3=psi
  	U=runif(k)
  	GRV=NULL
  	IND1=which(U<(p1))
  	IND2=which((U>(p1))&(U<(p1+p2)))
  	IND3=which(U>(p1+p2))
  	GRV[IND1]=-1
	GRV[IND2]=0
  	GRV[IND3]=1
  	return(GRV)
}


##################### RIS_RP ##########################

RIS_RP=function(x,y,x_test,alpha=alpha)
{
  a=0.02; b=0.02;
  z=standardize(x)
  p=ncol(x)
  n=nrow(x)
  yz=(y-(mean(y)*rep(1,n)))/sd(y)
  CR=apply(z,2,cor,y=yz)
  nu=2
  q_inn=(abs(CR)^nu)/max(abs(CR)^nu)  
  Mm=sample(seq(ceiling(2*log(p)),((3*n/4)-5)),105,replace=TRUE)
  Ppsi=runif(105,0.1,0.4)
  lgwght=NULL
  n.test=nrow(x_test)
  y_hat=lo.lim=up.lim=matrix(data=NA,nrow=105,ncol=n.test)
  for(i in 1:105)
  {
    m=Mm[i]
    psi=Ppsi[i]

# RIS part    
    Gmm=sapply(q_inn,rbinom,n=1,size=1) 
    I_use=which(Gmm==1)
    p_use=length(I_use)
    RM.star=matrix(data=rRM.MBCR(m*p_use,psi),nrow=p_use,ncol=m)

# Y_hat 
    S=t(x[,I_use]%*%RM.star)  	# mxn matrix of (XR)'
    SG=S%*%t(S)+diag(m)     	# mxm matrix of (I_m +(XR)'(XR))
    solve.SG=solve(SG)		# ((XR)'(XR)+I_m)^(-1)
    ZpZiZp=solve.SG%*%S     	# ((XR)'(XR)+I_m)^(-1)(XR)'
    beta_hat1=S%*%y		# (XR)'y
    Z.new=(x_test[,I_use]%*%RM.star)  # (X.new R)
    z.new.solve.SG=Z.new%*%solve.SG   #(X.new R)((XR)'(XR)+I_m)^(-1)
    y_hat[i,]=z.new.solve.SG%*%beta_hat1

# calculation of weight
    tS.solveSG.S=(diag(n)-t(S)%*%ZpZiZp)     # (I_n - XR ((XR)'(XR)+I_m)^(-1)(XR)')
    scale_par1=(t(y)%*%tS.solveSG.S%*%y)     # y'(I_n - XR ((XR)'(XR)+I_m)^(-1)(XR)')y

# paramterters and DF of t-distribution
    df=(n+2*a)
    const1=(scale_par1+2*b)/df

    for(kk in 1:n.test)
    {
		loc_par=y_hat[i,kk]
		scale_par=(1+t(Z.new[kk,])%*%solve.SG%*%(Z.new[kk,]))*const1
		pt=qt((0.5+alpha/2),lower.tail = TRUE,df=df,ncp=0)
		lo.lim[i,kk]=loc_par-(pt*sqrt(scale_par))
		up.lim[i,kk]=loc_par+(pt*sqrt(scale_par))
    }
  }
  Y.SA=colMeans(y_hat); lo.SA=colMeans(lo.lim); up.SA=colMeans(up.lim)
  Res.MBCR=cbind(Y.SA,lo.SA,up.SA)
  return(Res.MBCR)
}


##################### RIS-PCR #######################

RIS_PCR=function(x,y,x_test,alpha=alpha)
{
  a=0.02; b=0.02;
  z=standardize(x)
  p=ncol(x)
  n=nrow(x)
  yz=(y-(mean(y)*rep(1,n)))/sd(y)
  CR=apply(z,2,cor,y=yz)
  nu=2
  q_inn=(abs(CR)^nu)/max(abs(CR)^nu)  

  mse.min=1000000000
  Mm=sample(seq(ceiling(2*log(p)),((3*n/4)-5)),105,replace=TRUE)

  lgwght=NULL
  n.test=nrow(x_test)
  y_hat=lo.lim=up.lim=matrix(data=NA,nrow=105,ncol=n.test)
  for(i in 1:105)
  {
    m=Mm[i]
    Gmm=NULL
    Gmm=sapply(q_inn,rbinom,n=1,size=1) 
    I_use=which(Gmm==1)
    p_use=length(I_use)

    SVD=svd(x[,I_use])
    V=SVD$v
    D=SVD$d
    m=min(ncol(V),m)
    RM.star=V[,1:(min(ncol(V),m))] #

#    RM.star=matrix(data=rRM.MBCR(m*p,psi),nrow=p,ncol=m)

# Y_hat 
    m=ncol(RM.star)
    S=t(x[,I_use]%*%RM.star)  	# mxn matrix of (XR)'
    SG=S%*%t(S)+diag(m)     	# mxm matrix of (I_m +(XR)'(XR))
    solve.SG=solve(SG)		# ((XR)'(XR)+I_m)^(-1)
    ZpZiZp=solve.SG%*%S     	# ((XR)'(XR)+I_m)^(-1)(XR)'
    beta_hat1=S%*%y		# (XR)'y
    Z.new=(x_test[,I_use]%*%RM.star)  # (X.new R)
    z.new.solve.SG=Z.new%*%solve.SG   #(X.new R)((XR)'(XR)+I_m)^(-1)
    y_hat[i,]=z.new.solve.SG%*%beta_hat1

# calculation of weight
    tS.solveSG.S=(diag(n)-t(S)%*%ZpZiZp)     # (I_n - XR ((XR)'(XR)+I_m)^(-1)(XR)')
    scale_par1=(t(y)%*%tS.solveSG.S%*%y)     # y'(I_n - XR ((XR)'(XR)+I_m)^(-1)(XR)')y

# paramterters and DF of t-distribution
	df=(n+2*a)
	loc_par=y_hat[i,]
	const1=(scale_par1+2*b)/df
	scale_par=(diag(n.test)+z.new.solve.SG%*%t(Z.new))

	for(kk in 1:n.test)
	{
	  loc_par=y_hat[i,kk]
	  scale_par=(1+t(Z.new[kk,])%*%solve.SG%*%(Z.new[kk,]))*const1
	  pt=qt((0.5+alpha/2),lower.tail = TRUE,df=df,ncp=0)
	  lo.lim[i,kk]=loc_par-(pt*sqrt(scale_par))
	  up.lim[i,kk]=loc_par+(pt*sqrt(scale_par))
	}
	
  }
  Y.SA=colMeans(y_hat); lo.SA=colMeans(lo.lim); up.SA=colMeans(up.lim)
   Res.svd=cbind(Y.SA,lo.SA,up.SA)
  return(Res.svd)
}



