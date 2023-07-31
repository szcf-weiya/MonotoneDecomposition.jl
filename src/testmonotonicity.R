##
## download from https://www.stat.colostate.edu/~meyer/testmonotonicity.R
## 
# Code to do quadratic monotone spline test of mon vs not
#  Uses B-splines 
#
#  dir=0 if decreasing; otherwise will fit increasing function
#
#  matrix of covariates z does not contain one-vector and must
#  combine with the splines vectors to form linearly independent set
#
#  IF no covariates, set z=0
#
#  Scatterplot is (x[i],y[i]); k=number of knots (k-2 interior knots), 
#
# w is a vector of weights:
#           var(y[i]) is inversely proportional to w[i]

montest=function(x,y,z,w,k,nsim,dir){
	n=length(y)
	wmat=diag(sqrt(w))
	yw=wmat%*%y
	sm=1e-8
	spline=bspline(x,k)
	ns=length(spline$edges)/n
	if(length(z)==1){nz=0}else{nz=length(z)/n}
	m=ns+nz
	if(nz==1){z=matrix(z,ncol=1,nrow=n)}
	delta=matrix(nrow=m,ncol=n)
	if(nz>0){delta[1:nz,]=t(z)}
	delta[(nz+1):m,]=t(spline$edges)
	deltaw=delta%*%wmat
	qmat=deltaw%*%t(deltaw)
	umat=chol(qmat)
	smat=matrix(0,nrow=m,ncol=k)
	smat[(nz+1):m,]=spline$slopes
	if(dir==0){smat=-smat}
	xpl=spline$xpl
	bpl=spline$bpl
	uinv=solve(umat)
	amat=t(smat)%*%uinv
	rmat=matrix(runif(m*(m-k)),nrow=m-k)
	rmat=rmat-rmat%*%t(amat)%*%solve(amat%*%t(amat))%*%amat
	bmat=rbind(rmat,amat)	
	edges=t(solve(bmat))
	ysend=t(uinv)%*%deltaw%*%yw
	ans=conedualc(ysend,edges,nz)
	coef=ans$coef
	beta=uinv%*%t(edges)%*%coef
	yhat=t(delta)%*%beta
	qinv=solve(qmat)
	betau=qinv%*%deltaw%*%yw
	pmat=t(deltaw)%*%qinv%*%deltaw
	edf=sum(diag(pmat))
	ytil=t(delta)%*%betau
	ytilw=t(deltaw)%*%betau
	amat=t(smat)%*%qinv%*%deltaw
	minslope=min(amat%*%yw)
	tvec=spline$knots
	nk=length(tvec)
	if(minslope<sm){
		sighat=sum((yw-ytilw)^2)/(n-edf)
		vmat=t(smat)%*%qinv%*%smat*sighat
		vsqrt=chol(vmat)
		svec=t(smat)%*%beta
		slow=1:nsim
		for(isim in 1:nsim){
			zr=rnorm(nk)
			ssim=svec+t(vsqrt)%*%zr
#			ssim=t(vsqrt)%*%z
			slow[isim]=min(ssim)
		}
		pval=sum(slow<minslope)/nsim
	} else
	    pval=1
	
# find the pieces separately
	yxhat=t(delta[(nz+1):m,])%*%beta[(nz+1):m]
	yxtil=t(delta[(nz+1):m,])%*%betau[(nz+1):m]
	if(nz>1){
		yztil=t(delta[1:nz,])%*%betau[1:nz]
		yzhat=t(delta[1:nz,])%*%beta[1:nz]
	}else{yztil=delta[1,]*betau[1];yzhat=delta[1,]*beta[1]}
	ans=new.env()
	ans$pval=pval
	ans$cfit=yhat
	ans$ucfit=ytil
	ans$xpl=xpl
	ans$cpl=bpl%*%beta[(nz+1):m]
	ans$ucpl=bpl%*%(qinv%*%delta%*%y)[(nz+1):m]
	ans$beta=beta
	ans$betau=betau
	ans
}
	
##############################################################
# cone projection code: hinge algorithm in R
#  find theta-hat to minimize || y-theta ||^2
#  subject to theta = Xa + Db where b>=0
##############################################################

conedualc=function(y,delta,m0){
	n=length(y);m=length(delta)/n
	sm=1e-8;h=1:m<0;obs=1:m;check=0
	dm=1:m>0;dm[1:m0]=FALSE; h[1:m0]=TRUE
	if(m0==1){xmat=matrix(delta[1:m0,],nrow=n,ncol=1)}else{xmat=t(delta[1:m0,])}
	if (m0 == 0)
	    fit = 0
	else
    	fit=xmat%*%solve(t(xmat)%*%xmat)%*%t(xmat)%*%y
	b2=delta%*%(y-fit)
	if(max(b2)>sm){
		i=min(obs[b2==max(b2)])
		h[i]=TRUE
	}else{check=1;theta=1:n*0;a=0}
	while(check==0){
		xmat=matrix(delta[h,],ncol=n)
		a=solve(xmat%*%t(xmat))%*%xmat%*%y
		if(min(a)<(-sm)){
			avec=1:m*0;avec[h]=a
			i=min(obs[avec==min(avec)])
			h[i]=FALSE;check=0
		}else{
			check=1
			theta=t(xmat)%*%a
			b2=delta%*%(y-theta)
			if(max(b2)>sm){
				i=min(obs[b2==max(b2)])		
				h[i]=TRUE;check=0
			}
		}
	}
	ans=new.env()
	ans$fit=theta
	b=1:m*0;b[h]=a
	ans$coef=b
	ans$df=sum(h)
	ans
}

########################################
#       MAKE THE EDGE VECTORS          #
########################################

##############################################################
# B-spline quadratic basis
# returns basis functions and slopes of basis functions
# at the knots
bspline=function(x,m){
	tk=0:(m-1)/(m-1)*(max(x)-min(x))+min(x)
	k=3
	t=1:(m+2*k-2)*0
	t[1:(k-1)]=min(x);t[(m+k):(m+2*k-2)]=max(x)
	t[k:(m+k-1)]=tk
	n=length(x)
	sm=1e-8
	h=t[4]-t[3]

	bmat=matrix(1:(n*(m+k-2))*0,nrow=n)
	index=x>=t[3]&x<=t[4]
	bmat[index,1]=(t[4]-x[index])^2
	bmat[index,2]=2*(x[index]-t[2])*(t[4]-x[index])+(t[5]-x[index])*(x[index]-t[3])
	index=x>=t[4]&x<=t[5]
	bmat[index,2]=(t[5]-x[index])^2
	for( j in 3:(m-1) ){
		index=x>=t[j]&x<=t[j+1]
		bmat[index,j]=(x[index]-t[j])^2
		index=x>=t[j+1]&x<=t[j+2]
		bmat[index,j]=(x[index]-t[j])*(t[j+2]-x[index])+(x[index]-t[j+1])*(t[j+3]-x[index])
		index=x>=t[j+2]&x<=t[j+3]
		bmat[index,j]=(t[j+3]-x[index])^2
	}
	index=x>=t[m]&x<=t[m+1]
	bmat[index,m]=(x[index]-t[m])^2
	index=x>=t[m+1]&x<=t[m+2]
	bmat[index,m]=(x[index]-t[m])*(t[m+2]-x[index])+2*(x[index]-t[m+1])*(t[m+3]-x[index])
	index=x>=t[m+1]&x<=t[m+2]
	bmat[index,m+1]=(x[index]-t[m+1])^2

# plotting splines

	xpl=0:1000/1000*(max(x)-min(x))+min(x)
	bpl=matrix(1:(1001*(m+k-2))*0,nrow=1001)
	index=xpl>=t[3]&xpl<=t[4]
	bpl[index,1]=(t[4]-xpl[index])^2
	bpl[index,2]=2*(xpl[index]-t[2])*(t[4]-xpl[index])+(t[5]-xpl[index])*(xpl[index]-t[3])
	index=xpl>=t[4]&xpl<=t[5]
	bpl[index,2]=(t[5]-xpl[index])^2
	for( j in 3:(m-1) ){
		index=xpl>=t[j]&xpl<=t[j+1]
		bpl[index,j]=(xpl[index]-t[j])^2
		index=xpl>=t[j+1]&xpl<=t[j+2]
		bpl[index,j]=(xpl[index]-t[j])*(t[j+2]-xpl[index])+(xpl[index]-t[j+1])*(t[j+3]-xpl[index])
		index=xpl>=t[j+2]&xpl<=t[j+3]
		bpl[index,j]=(t[j+3]-xpl[index])^2
	}
	index=xpl>=t[m]&xpl<=t[m+1]
	bpl[index,m]=(xpl[index]-t[m])^2
	index=xpl>=t[m+1]&xpl<=t[m+2]
	bpl[index,m]=(xpl[index]-t[m])*(t[m+2]-xpl[index])+2*(xpl[index]-t[m+1])*(t[m+3]-xpl[index])
	index=xpl>=t[m+1]&xpl<=t[m+2]
	bpl[index,m+1]=(xpl[index]-t[m+1])^2

# matrix of slopes

	slopes=matrix(0,ncol=m,nrow=m+k-2)
	slopes[1,1]=-2*h
	slopes[m+k-2,m]=2*h
	slopes[2,1]=4*h
	slopes[2,2]=-2*h
	slopes[m+k-3,m]=-4*h
	slopes[m+k-3,m-1]=2*h
	if(m>=4){
		for(j in 3:(m+k-4)){
			slopes[j,j-1]=2*h
			slopes[j,j]=-2*h
		}
	}

	bmat[,1]=bmat[,1]*2
	bmat[,m+1]=bmat[,m+1]*2
	slopes[1,]=slopes[1,]*2
	slopes[m+1,]=slopes[m+1,]*2
	bpl[,1]=bpl[,1]*2
	bpl[,m+1]=bpl[,m+1]*2
	mb=max(bpl)
	slopes=slopes/mb
	bpl=bpl/mb
	bmat=bmat/mb
	ans=new.env()
	ans$edges=bmat
	ans$slopes=slopes
	ans$knots=tk
	ans$xpl=xpl
	ans$bpl=bpl
	ans
}

