## Method 1 : ISIS+SCAD
# ISIS.SCAD=function(x,y,x_test,y_test)
# {
#   n.test=nrow(x_test)
#   n=nrow(x)
#   nt=nrow(x_test)
# 	FIT=SIS(x,y,family="gaussian",penalty="SCAD",iter=TRUE) 
# 	S_isis=FIT$ix
#   if(length(S_isis)==0) {S_isis=FIT$ix0}
# 	x_select=x[,S_isis] 
#   new.dat=x_test[,S_isis]
# rownames(x_select)=NULL
# rownames(new.dat)=NULL
# olddata=data.frame(c(y,y_test),rbind(x_select,new.dat))
# regress=olddata[1:n,(2:ncol(olddata))]
# FIT2=lm(olddata[1:n,1]~as.matrix(regress))
# regress = olddata[((n+1):(n+nt)),(2:ncol(olddata))]
# Pred=predict(FIT2,regress,interval="prediction",level=0.5)
#   Pred1=predict(FIT2,as.data.frame(regress),type="response")
#   mse.scad=mean((as.vector(y_test)-as.vector(Pred1))^2)
# 
# lo=Pred[,2]
# up=Pred[,3]
# dst=up-lo
#   lower=sign(y_test-lo)
#   upper=sign(up-y_test)
#   
#   count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
# return(c(mse.scad,count,mean(dst)))
# }
# 
# ## Method 2 : ISIS+MCP  
# ISIS.MCP=function(x,y,x_test,y_test)
# {
#   n.test=nrow(x_test)
# n=nrow(x)
# nt=nrow(x_test)
# 	FIT=SIS(x,y,family="gaussian",penalty="MCP",iter=TRUE) 
# 	S_isis=FIT$ix
# if(length(S_isis)==0) {S_isis=FIT$ix0}
# 	x_select=x[,S_isis]    
# 
# new.dat=x_test[,S_isis]
# rownames(x_select)=NULL
# rownames(new.dat)=NULL
# 
#   olddata=data.frame(c(y,y_test),rbind(x_select,new.dat))
#   regress=olddata[1:n,(2:ncol(olddata))]
#   FIT2=lm(olddata[1:n,1]~as.matrix(regress))
#   regress = olddata[((n+1):(n+nt)),(2:ncol(olddata))]
#   Pred=predict(FIT2,regress,interval="prediction",level=0.5)
#   Pred1=predict(FIT2,as.data.frame(regress),type="response")
#   mse.mcp=mean((as.vector(y_test)-as.vector(Pred1))^2)
# 
#   lo=Pred[,2]
#   up=Pred[,3]
#   dst=(up-lo)
#   lower=sign(y_test-lo)
#   upper=sign(up-y_test)
#   count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
# return(c(mse.mcp,count,mean(dst)))
# }
# 

## Method 3 : BigLASSO, screening : sequential strong rule (SSR), "Strong rules for discarding predictors in lasso-type problems, 2011, JRSSB"
BigLasso=function(x,y,x_test,y_test)
{
   n.test=nrow(x_test)
   n=nrow(x)
   nt=nrow(x_test)

	X.bm <- as.big.matrix(x)
#	L=seq(0.0005,2,length=200)
	fit=biglasso(X.bm,y,family ='gaussian',penalty = "lasso",ncores=1,alpha=1)
	L=fit$lambda
	k=apply(fit$beta,2,function(u){length(which(u!=0))})
	BIC=n*log(fit$loss/n)+k*log(n)
  AIC=n*log(fit$loss/n)+2*k
	l.opt=fit$lambda[which.min(AIC)]
	print(l.opt)

	fit=glmnet(x,y,family ='gaussian',alpha=1,lambda=l.opt,standardize = FALSE)
	y.fit=predict(fit,x_test, type="response")	
	mse.lasso=mean((as.vector(y.fit)-y_test)^2)
	
	u=NULL
	for(j in 1:n)
	{
	  y.train=y[-j]; y.test=y[j];
	  x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
	  lasso.mod <- glmnet(x.train, y.train, alpha = 1, lambda=l.opt,family="gaussian", intercept=FALSE,standardize = FALSE)
	  y.hat <- predict(lasso.mod, newx = x.test,type="response")
	  u[j]=y.test-y.hat
	}
	u1=quantile(u,0.25)
	u2=quantile(u,0.75)
	lo=y.fit+u1
	up=y.fit+u2
	dst=(up-lo)
	lower=sign(y_test-as.numeric(lo))
	upper=sign(as.numeric(up)-y_test)
	count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
	c(mse.lasso,count,mean(dst))
	return(c(mse.lasso,count,mean(dst)))
	
}

## Method 3 : CV BigLASSO, screening : sequential strong rule (SSR), "Strong rules for discarding predictors in lasso-type problems, 2011, JRSSB"
ISIS.SCAD=function(x,y,x_test,y_test)   ## CAREFUL!!
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
  
  #X.bm <- as.big.matrix(x)
  fit=cv.glmnet(x,y,family ='gaussian')
#  k=apply(fit$beta,2,function(u){length(which(u!=0))})
#  BIC=n*log(fit$loss/n)+k*log(n)
#  AIC=n*log(fit$loss/n)+2*k
#  l.opt=fit$lambda[which.min(AIC)]
#  print(l.opt)
  
#  fit=glmnet(x,y,family ='gaussian',alpha=1,lambda=l.opt,standardize = FALSE)
  y.fit=predict(fit,x_test, type="response")	
  mse.lasso=mean((as.vector(y.fit)-y_test)^2)
  
  u=NULL
  for(j in 1:n)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
    lasso.mod <- cv.glmnet(x.train, y.train, family="gaussian", intercept=FALSE,standardize = FALSE)
    y.hat <- predict(lasso.mod, newx = x.test,type="response")
    u[j]=y.test-y.hat
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  upper=sign(as.numeric(up)-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  c(mse.lasso,count,mean(dst))
  return(c(mse.lasso,count,mean(dst)))
  
}


## Method 4 : Ridge
BigRidge=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
  X.bm <- as.big.matrix(x)
  X.bm <- as.big.matrix(x)
  #L=seq(0.001,10,length=250)
#  L=seq(0.0005,2,length=200)
  fit=biglasso(X.bm,y,family ='gaussian',ncores=1,alpha=0)
  L=fit$lambda
  k=apply(fit$beta,2,function(u){length(which(u!=0))})
  #BIC=n*log(fit$loss/n)+k*log(n)
  AIC=n*log(fit$loss/n)+2*k
  l.opt=fit$lambda[which.min(AIC)]
  print(l.opt)
  fit=glmnet(x,y,family ='gaussian',alpha=0,lambda=l.opt,standardize = FALSE)
  y.fit=predict(fit,x_test, type="response")	
  mse.lasso=mean((as.vector(y.fit)-y_test)^2)
  u=NULL
  for(j in 1:n)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
    lasso.mod <- glmnet(x.train, y.train, alpha = 0, lambda=l.opt,family="gaussian", intercept=FALSE,standardize = FALSE)
    y.hat <- predict(lasso.mod, newx = x.test,type="response")
    u[j]=y.test-y.hat
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  upper=sign(as.numeric(up)-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  c(mse.lasso,count,mean(dst))
  return(c(mse.lasso,count,mean(dst)))

}


## Method 4 : Elastic net
BigENet=function(x,y,x_test,y_test)
{
 	n.test=nrow(x_test)
	n=nrow(x)
	nt=nrow(x_test)
	X.bm <- as.big.matrix(x)
	X.bm <- as.big.matrix(x)
	# L=seq(0.001,10,length=250)
	L=seq(0.0005,2,length=200)
#	fit=biglasso(X.bm,y,family ='gaussian',penalty ="enet",ncores=1,alpha=0.5,lambda=L)
	fit=biglasso(X.bm,y,family ='gaussian',penalty = "enet",ncores=1,alpha=0.5)
  L=fit$lambda
	k=apply(fit$beta,2,function(u){length(which(u!=0))})
#	BIC=n*log(fit$loss/n)+k*log(n)
AIC=n*log(fit$loss/n)+2*k
	l.opt=fit$lambda[which.min(AIC)]
  print(l.opt)
  fit=glmnet(x,y,family ='gaussian',alpha=0.5,lambda=l.opt,standardize = FALSE)
  y.fit=predict(fit,x_test, type="response")	
  mse.lasso=mean((as.vector(y.fit)-y_test)^2)
  
  u=NULL
  for(j in 1:n)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
    lasso.mod <- glmnet(x.train, y.train, alpha = 0.5, lambda=l.opt,family="gaussian", intercept=FALSE,standardize = FALSE)
    y.hat <- predict(lasso.mod, newx = x.test,type="response")
    u[j]=y.test-y.hat
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  	upper=sign(as.numeric(up)-y_test)
  	count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  	c(mse.lasso,count,mean(dst))
	  
  	return(c(mse.lasso,count,mean(dst)))
}



## Method 7 : SVD
SVD=function(x,y,x_test,y_test)
{
 n.test=nrow(x_test)
n=nrow(x)
nt=nrow(x_test)
#x=standardize(x)
#x_test=standardize(x_test)

	p=ncol(x)
	SVD=fast.svd(x)
#	SVD=ssvd(x)
	D=SVD$d
	RM.star=SVD$v
	col.rm=ncol(RM.star)
	if(col.rm>20)
	{ RM.star=RM.star[,(1:(col.rm-2))] }
	
#	beta_hat=solve(t(x%*%RM.star)%*%(x%*%RM.star))%*%t(x%*%RM.star)%*%y
# y.fit=(x_test%*%RM.star)%*%beta_hat
#	mse.RM=mean((as.vector(y_test)-y.fit)^2)

	x_select=x%*%RM.star
	new.dat=x_test%*%RM.star
	olddata=data.frame(c(y,y_test),rbind(x_select,new.dat))
	regress=olddata[1:n,(2:ncol(olddata))]
	FIT2=lm(olddata[1:n,1]~as.matrix(regress))
	regress = olddata[((n+1):(n+n.test)),(2:ncol(olddata))]
	Pred1=predict(FIT2,as.data.frame(regress),type="response")
	Pred=predict(FIT2,as.data.frame(regress),interval="prediction",level=0.5)
	lo=Pred[,2]
	up=Pred[,3]
	dst=up-lo
	lower=sign(y_test-lo)
	upper=sign(up-y_test)
	mse.RM=mean((as.vector(y_test)-as.vector(Pred1))^2)
	count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
	return(c(mse.RM,count,mean(dst)))
}
# 	df=nrow(x)
# 	y.sd=(sd(y.fit))^2
# 	scale.par=sqrt(y.sd+mse.RM)
# 	pt=qt(0.75,lower.tail = TRUE,df=df,ncp=0)
# 	lo=y.fit-(pt*scale.par)
# 	up=y.fit+(pt*scale.par)
# dst=(up-lo)
#   lower=sign(y_test-as.numeric(lo))
#   upper=sign(as.numeric(up)-y_test)
#   count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
# 	return(c(count,mean(dst)))


# 
# ### CV lasso
# BigCVLASSO=function(x,y,x_test,y_test)
# {
#   t1=Sys.time()
#   X.bm <- as.big.matrix(x)
#   fit=cv.biglasso(X.bm,t(y),family ='gaussian',ncores=1)
#   x.t=as.big.matrix(x_test)
#   y.fit=predict(fit,x.t, type="response")	
#   mse.cvlas=mean((as.vector(y.fit)-y_test)^2)
#   t2=difftime(Sys.time(),t1,units="mins")
#   Res.cvlas=c(mse.cvlas,t2)
#   return(Res.cvlas)
# }
# 


## Method 8 : Sparse PCR

SSVD=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
 # x=standardize(x)
 # x_test=standardize(x_test)
  col.rm=20
  p=ncol(x)
  SVD=SPC(x,K=col.rm)
  RM.star=SVD$v
  col.rm=ncol(RM.star)
  if(col.rm>20)
  { RM.star=RM.star[,(1:(col.rm-2))] }

  x_select=x%*%RM.star
  new.dat=x_test%*%RM.star
  olddata=data.frame(c(y,y_test),rbind(x_select,new.dat))
  regress=olddata[1:n,(2:ncol(olddata))]
  FIT2=lm(olddata[1:n,1]~as.matrix(regress))
  regress = olddata[((n+1):(n+n.test)),(2:ncol(olddata))]
  Pred1=predict(FIT2,as.data.frame(regress),type="response")
  Pred=predict(FIT2,as.data.frame(regress),interval="prediction",level=0.5)
  lo=Pred[,2]
  up=Pred[,3]
  dst=up-lo
  lower=sign(y_test-lo)
  upper=sign(up-y_test)
  mse.RM=mean((as.vector(y_test)-as.vector(Pred1))^2)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  return(c(mse.RM,count,mean(dst)))
  
}

## Method 8 : RPCA

RPCA=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
 # x=standardize(x)
 # x_test=standardize(x_test)
  col.rm=20
  p=ncol(x)
  SVD=rsvd(x,k=col.rm)
  RM.star=SVD$v
  col.rm=ncol(RM.star)
  if(col.rm>20)
  { RM.star=RM.star[,(1:(col.rm-2))] }

  x_select=x%*%RM.star
  new.dat=x_test%*%RM.star
  olddata=data.frame(c(y,y_test),rbind(x_select,new.dat))
  regress=olddata[1:n,(2:ncol(olddata))]
  FIT2=lm(olddata[1:n,1]~as.matrix(regress))
  regress = olddata[((n+1):(n+n.test)),(2:ncol(olddata))]
  Pred1=predict(FIT2,as.data.frame(regress),type="response")
  Pred=predict(FIT2,as.data.frame(regress),interval="prediction",level=0.5)
  lo=Pred[,2]
  up=Pred[,3]
  dst=up-lo
  lower=sign(y_test-lo)
  upper=sign(up-y_test)
  mse.RM=mean((as.vector(y_test)-as.vector(Pred1))^2)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  return(c(mse.RM,count,mean(dst)))
  
}

BASAD=function(x,y,x_test,y_test)
{
  nt=n.test=nrow(x_test)
  n=nrow(x)
  # x=standardize(x)
  # x_test=standardize(x_test)
  #fit=basad(x=x,y=y,select.cri="BIC",BIC.maxsize=50)
fit=basad(x=x,y=y,select.cri="AIC",AIC.maxsize=50)
  Pred1=predict(fit,testx=x_test)
  mse.basad=mean((as.vector(y_test)-as.vector(Pred1))^2)

  # Prediction
  S_isis=which(fit$modelZ!=0)
  x_select=cbind(rep(1,n),x)
  if(length(intersect(as.numeric(S_isis),(1:ncol(x_select))))!=length(S_isis))
  {
    uv=which((S_isis%in%(1:ncol(x_select)))==TRUE)
    S_isis=S_isis[uv]
  }
  x_select=x[,as.numeric(S_isis)]
  new.dat=cbind(rep(1,nt),x_test)
  new.dat=x_test[,as.numeric(S_isis)]
  rownames(x_select)=NULL
  rownames(new.dat)=NULL
  olddata=data.frame(c(y,y_test),rbind(x_select,new.dat))
  regress=olddata[1:n,(2:ncol(olddata))]
  FIT2=lm(olddata[1:n,1]~as.matrix(regress))
  regress = olddata[((n+1):(n+nt)),(2:ncol(olddata))]
  Pred=predict(FIT2,regress,interval="prediction",level=0.5)
  lo=Pred[,2]
  up=Pred[,3]
  dst=up-lo
  lower=sign(y_test-lo)
  upper=sign(up-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  return(c(mse.basad,count,mean(dst)))
}

## SSLASSo :
SSlasso=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
# usual 
#  fit=SSLASSO::SSLASSO(x,y,lambda1=1,lambda0=seq(1,100,1),penalty = "adaptive"); 
# for eye data only
  fit=SSLASSO::SSLASSO(x,y,lambda1=1,lambda0=(1+(seq(1,100,1)/20)),penalty = "adaptive"); 
  S_isis=fit$model
  x_select=x[,S_isis] 
  new.dat=x_test[,S_isis]
  rownames(x_select)=NULL
  rownames(new.dat)=NULL
  olddata=data.frame(c(y,y_test),rbind(x_select,new.dat))
  regress=olddata[1:n,(2:ncol(olddata))]
  FIT2=lm(olddata[1:n,1]~as.matrix(regress))
  regress = olddata[((n+1):(n+nt)),(2:ncol(olddata))]
  Pred=predict(FIT2,regress,interval="prediction",level=0.5)
  Pred1=predict(FIT2,as.data.frame(regress),type="response")
  mse.scad=mean((as.vector(y_test)-as.vector(Pred1))^2)
  
  lo=Pred[,2]
  up=Pred[,3]
  dst=up-lo
  lower=sign(y_test-lo)
  upper=sign(up-y_test)
  
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  return(c(mse.scad,count,mean(dst)))
}


# Method 0: One-Step SCAD
OneStepScad=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
  
  scad.out<-linreg.illa(x=x, y=as.matrix(y), xtune=x, ytune=as.matrix(y), penalty="scad", nlambda=200, lambda.min.ratio=0.0005)

  l.opt=scad.out$lambda.opt
  intercept=scad.out$coef.icpt
  bbeta=scad.out$coef.est
  S_isis=which(bbeta!=0)
#	x_select=cbind(rep(1,n),x[,S_isis])
  new.dat=cbind(rep(1,n.test),x_test[,S_isis])
#print(dim(new.dat))
#print(length(bbeta[S_isis]))
  y.fit=Pred1=new.dat%*%c(intercept,bbeta[S_isis])
  mse.scad=mean((as.vector(y_test)-as.vector(Pred1))^2)
  u=NULL
  J=sample(1:n,min(50,nrow(x)),replace=FALSE); i=1
  for(j in J)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=as.vector(x[j,]);
    scad.out<-linreg.illa(x=x.train, y=as.matrix(y.train), xtune=x.train, ytune=as.matrix(y.train), penalty="scad", nlambda=2, lambda.min.ratio=l.opt)
    intercept=scad.out$coef.icpt
    bbeta=scad.out$coef.est
    S_isis=which(bbeta!=0)
#    x_select=cbind(rep(1,n),x[,S_isis])
    new.dat=t(as.matrix(c(1,x.test[S_isis])))
#    print(dim(new.dat))
#    print(length(bbeta[S_isis]))
    
    y.hat=new.dat%*%c(intercept,bbeta[S_isis])
    u[i]=y.test-y.hat
    i=i+1
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  upper=sign(as.numeric(up)-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  return(c(mse.scad,count,mean(dst)))
}


################################ MCP / SCAD (ncpen) ###############################


SCAD.ncpen.L=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
  library(ncpen)
  L=seq(0.0005,1,length=200)
  FIT=ncpen(y.vec=y,x.mat=x,family="gaussian",penalty="scad",lambda=L)
  Choose_lambda=gic.ncpen(FIT, weight =2)
  l.opt=Choose_lambda$opt.lambda
  FIT=ncpen(y.vec=y,x.mat=x,family="gaussian",penalty="scad",lambda=l.opt)
  y.fit.mat=predict(FIT,type="y",new.x.mat=x_test)
  opt.y=which(FIT$lambda==l.opt)
  if(length(opt.y)==0){opt.y=length(FIT$lambda)}
  y.fit=y.fit.mat[,opt.y]
  mse.lasso=mean((as.vector(y.fit)-y_test)^2)

  u=NULL
  for(j in 1:n)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
    lasso.mod <- ncpen(y.train, x.train, penalty="scad",family="gaussian")
    Choose_lambda=gic.ncpen(lasso.mod, weight =2,verbose = FALSE)
    l.opt=Choose_lambda$opt.lambda
    y.hat.mat <- predict(lasso.mod, new.x.mat = x.test,type="y")
    opt.y=which(lasso.mod$lambda==l.opt)
    if(length(opt.y)==0){opt.y=length(lasso.mod$lambda)}
    y.hat=y.hat.mat[,opt.y]
    u[j]=y.test-y.hat
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  upper=sign(as.numeric(up)-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  c(mse.lasso,count,mean(dst))
  return(c(mse.lasso,count,mean(dst)))
}


SCAD.ncpen.A=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
  library(ncpen)
  #L=seq(0.0005,1,length=200)
  FIT=ncpen(y.vec=y,x.mat=x,family="gaussian",penalty="scad")
  Choose_lambda=gic.ncpen(FIT, weight =2)
  l.opt=Choose_lambda$opt.lambda
  FIT=ncpen(y.vec=y,x.mat=x,family="gaussian",penalty="scad",lambda=l.opt)
  y.fit.mat=predict(FIT,type="y",new.x.mat=x_test)
  opt.y=which(FIT$lambda==l.opt)
  if(length(opt.y)==0){opt.y=length(FIT$lambda)}
  y.fit=y.fit.mat[,opt.y]
  mse.lasso=mean((as.vector(y.fit)-y_test)^2)
  
  u=NULL
  for(j in 1:n)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
    lasso.mod <- ncpen(y.train, x.train, penalty="scad",family="gaussian")
    Choose_lambda=gic.ncpen(lasso.mod, weight =2,verbose = FALSE)
    l.opt=Choose_lambda$opt.lambda
    y.hat.mat <- predict(lasso.mod, new.x.mat = x.test,type="y")
    opt.y=which(lasso.mod$lambda==l.opt)
    if(length(opt.y)==0){opt.y=length(lasso.mod$lambda)}
    y.hat=y.hat.mat[,opt.y]
    u[j]=y.test-y.hat
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  upper=sign(as.numeric(up)-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  c(mse.lasso,count,mean(dst))
  return(c(mse.lasso,count,mean(dst)))
}

MCP.ncpen.L=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
  library(ncpen)
  L=seq(0.0005,1,length=200)
  FIT=ncpen(y.vec=y,x.mat=x,family="gaussian",penalty="mcp",lambda=L)
  Choose_lambda=gic.ncpen(FIT, weight =2,verbose = FALSE)
  l.opt=Choose_lambda$opt.lambda
  FIT=ncpen(y.vec=y,x.mat=x,family="gaussian",penalty="mcp",lambda=l.opt)
  y.fit.mat=predict(FIT,type="y",new.x.mat=x_test)
  opt.y=which(FIT$lambda==l.opt)
  if(length(opt.y)==0){opt.y=length(FIT$lambda)}
  y.fit=y.fit.mat[,opt.y]
  mse.lasso=mean((as.vector(y.fit)-y_test)^2)

  u=NULL
  for(j in 1:n)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
    lasso.mod <- ncpen(y.train, x.train, penalty="mcp",family="gaussian")
    Choose_lambda=gic.ncpen(lasso.mod, weight =2,verbose = FALSE)
    l.opt=Choose_lambda$opt.lambda
    y.hat.mat <- predict(lasso.mod, new.x.mat = x.test,type="y")
    opt.y=which(lasso.mod$lambda==l.opt)
    if(length(opt.y)==0){opt.y=length(lasso.mod$lambda)}
    y.hat=y.hat.mat[,opt.y]
    u[j]=y.test-y.hat
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  upper=sign(as.numeric(up)-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  c(mse.lasso,count,mean(dst))
  return(c(mse.lasso,count,mean(dst)))
}


MCP.ncpen.A=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
  library(ncpen)
  L=seq(0.0005,1,length=200)
  FIT=ncpen(y.vec=y,x.mat=x,family="gaussian",penalty="mcp")
  Choose_lambda=gic.ncpen(FIT, weight =2,verbose = FALSE)
  l.opt=Choose_lambda$opt.lambda
  FIT=ncpen(y.vec=y,x.mat=x,family="gaussian",penalty="mcp",lambda=l.opt)
  y.fit.mat=predict(FIT,type="y",new.x.mat=x_test)
  opt.y=which(FIT$lambda==l.opt)
  if(length(opt.y)==0){opt.y=length(FIT$lambda)}
  y.fit=y.fit.mat[,opt.y]
  mse.lasso=mean((as.vector(y.fit)-y_test)^2)
  
  u=NULL
  for(j in 1:n)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
    lasso.mod <- ncpen(y.train, x.train, penalty="mcp",family="gaussian")
    Choose_lambda=gic.ncpen(lasso.mod, weight =2,verbose = FALSE)
    l.opt=Choose_lambda$opt.lambda
    y.hat.mat <- predict(lasso.mod, new.x.mat = x.test,type="y")
    opt.y=which(lasso.mod$lambda==l.opt)
    if(length(opt.y)==0){opt.y=length(lasso.mod$lambda)}
    y.hat=y.hat.mat[,opt.y]
    u[j]=y.test-y.hat
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  upper=sign(as.numeric(up)-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  c(mse.lasso,count,mean(dst))
  return(c(mse.lasso,count,mean(dst)))
}


################################ MCP / SCAD (ncpen) ###############################

SCAD.ncvreg.L=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
  library(ncvreg)
  L=seq(0.0005,1,length=200)
  FIT=ncvreg(X=x,y=y,family="gaussian",penalty="SCAD",lambda = L)
	k=apply(FIT$beta,2,function(u){length(which(u!=0))})
  AIC=n*log(FIT$loss/n)+2*k
	l.opt=FIT$lambda[which.min(AIC)]

  FIT=ncvreg(X=x,y=y,family="gaussian",penalty="SCAD",lambda=l.opt)
  y.fit=predict(FIT,type="response",X=x_test)	
  mse.lasso=mean((as.vector(y.fit)-y_test)^2)
  
  u=NULL
  for(j in 1:n)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
    lasso.mod <- ncvreg( x.train,y.train, penalty="SCAD",family="gaussian",lambda=l.opt)
    y.hat <- predict(lasso.mod, X = x.test,type="response")
    u[j]=y.test-y.hat
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  upper=sign(as.numeric(up)-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  c(mse.lasso,count,mean(dst))
  return(c(mse.lasso,count,mean(dst)))
}



SCAD.ncvreg.A=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
  library(ncvreg)
  #L=seq(0.0005,1,length=200)
  FIT=ncvreg(X=x,y=y,family="gaussian",penalty="SCAD")
  k=apply(FIT$beta,2,function(u){length(which(u!=0))})
  AIC=n*log(FIT$loss/n)+2*k
  l.opt=FIT$lambda[which.min(AIC)]
  
  FIT=ncvreg(X=x,y=y,family="gaussian",penalty="SCAD",lambda=l.opt)
  y.fit=predict(FIT,type="response",X=x_test)	
  mse.lasso=mean((as.vector(y.fit)-y_test)^2)
  
  u=NULL
  for(j in 1:n)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
    lasso.mod <- ncvreg( x.train,y.train, penalty="SCAD",family="gaussian",lambda=l.opt)
    y.hat <- predict(lasso.mod, X = x.test,type="response")
    u[j]=y.test-y.hat
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  upper=sign(as.numeric(up)-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  c(mse.lasso,count,mean(dst))
  return(c(mse.lasso,count,mean(dst)))
}



MCP.ncvreg.L=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
  library(ncvreg)
  L=seq(0.0005,1,length=200)
  FIT=ncvreg(X=x,y=y,family="gaussian",penalty="MCP",lambda=L)
  k=apply(FIT$beta,2,function(u){length(which(u!=0))})
  AIC=n*log(FIT$loss/n)+2*k
  l.opt=FIT$lambda[which.min(AIC)]
  
  FIT=ncvreg(X=x,y=y,family="gaussian",penalty="MCP",lambda=l.opt)
  y.fit=predict(FIT,type="response",X=x_test)	
  mse.lasso=mean((as.vector(y.fit)-y_test)^2)
  
  u=NULL
  for(j in 1:n)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
    lasso.mod <- ncvreg( x.train,y.train, penalty="MCP",family="gaussian",lambda=l.opt)
    y.hat <- predict(lasso.mod, X = x.test,type="response")
    u[j]=y.test-y.hat
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  upper=sign(as.numeric(up)-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  c(mse.lasso,count,mean(dst))
  return(c(mse.lasso,count,mean(dst)))
}


MCP.ncvreg.A=function(x,y,x_test,y_test)
{
  n.test=nrow(x_test)
  n=nrow(x)
  nt=nrow(x_test)
  library(ncvreg)
  #L=seq(0.0005,1,length=200)
  FIT=ncvreg(X=x,y=y,family="gaussian",penalty="MCP")
  k=apply(FIT$beta,2,function(u){length(which(u!=0))})
  AIC=n*log(FIT$loss/n)+2*k
  l.opt=FIT$lambda[which.min(AIC)]
  
  FIT=ncvreg(X=x,y=y,family="gaussian",penalty="MCP",lambda=l.opt)
  y.fit=predict(FIT,type="response",X=x_test)	
  mse.lasso=mean((as.vector(y.fit)-y_test)^2)
  
  u=NULL
  for(j in 1:n)
  {
    y.train=y[-j]; y.test=y[j];
    x.train=x[-j,]; x.test=t(as.matrix(x[j,]));
    lasso.mod <- ncvreg( x.train,y.train, penalty="MCP",family="gaussian",lambda=l.opt)
    y.hat <- predict(lasso.mod, X = x.test,type="response")
    u[j]=y.test-y.hat
  }
  u1=quantile(u,0.25)
  u2=quantile(u,0.75)
  lo=y.fit+u1
  up=y.fit+u2
  dst=(up-lo)
  lower=sign(y_test-as.numeric(lo))
  upper=sign(as.numeric(up)-y_test)
  count=(n.test-(length(which(lower==-1))+length(which(upper==-1))))/n.test
  c(mse.lasso,count,mean(dst))
  return(c(mse.lasso,count,mean(dst)))
}
