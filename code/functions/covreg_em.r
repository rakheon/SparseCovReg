# The functions for implementing Hoff and Niu (2012).
ldmvnorm<-function(y,mu=rep(0,length(y)),Sig=diag(1,length(y)))
{
  -.5*( length(y)*log(2*pi) + log(det(Sig)) + t(y-mu)%*%solve(Sig)%*%(y-mu)  )
}
covreg.em<-function(Y,X,R=dim(Y)[2],tol=1e-10,itmax=1000)
{

  p<-dim(Y)[2] ; q<-dim(X)[2] ; n<-dim(Y)[1]
  S0<-cov(Y);iS0<-solve(S0);B<-matrix(rnorm(p*q*R),p,q*R)*1e-4;iter<-0
  LL<-NULL
  rll<-10
  while( rll > tol  & iter<itmax)
  {
    B0<-B ; iter<-iter+1

    ### find expectation, var of z
    Vz<-array(dim=c(R,R,n)) ; Mz<-matrix(nrow=n,ncol=R)
    for(i in 1:n)
    {
      Bx<-apply(array(B,dim=c(p,q,R)),3,"%*%",X[i,])
      Vz[,,i]<-solve(  t(Bx)%*%iS0%*%Bx + diag(R) )
      Mz[i,]<-Vz[,,i]%*%t(Bx)%*%iS0%*%Y[i,]
    }
    ###

    ### obtain MLEs
    Y1<-Y ; X1<-NULL ; for(r in 1:R) { X1<-cbind(X1,diag(Mz[,r])%*%X  )}
    Y0<-matrix(0,nrow=n*R,ncol=p) ; X0<-NULL
    for(i in 1:n)
    {
      xi<-matrix(outer(X[i,],diag(R)),nrow=R*q,ncol=R)
      ZZ<-xi%*%Vz[,,i]%*%t(xi) ; ZZ<-.5*(ZZ+t(ZZ))
      Z<-eigen(ZZ);Z<-Z$vec[,1:R]%*%diag(sqrt(Z$val[1:R]),nrow=R) 
      X0<-rbind(X0,t(Z))
    }
    YA<-rbind(Y0,Y1) ; XA<-rbind(X0,X1)

    B<-t(YA)%*%XA%*%solve(t(XA)%*%XA) 
    E<-YA-XA%*%t(B)
    S0<- (t(E)%*%E)/dim(Y)[1]  ; iS0<-solve(S0)
    ###


    ###
    if(iter%%5==0)
    {
      ll<-0
      for(i in 1:dim(Y)[1])
      {
        xi<-matrix(outer(X[i,],diag(R)),nrow=R*q,ncol=R)
        ll<-ll+ldmvnorm(Y[i,],Sig=S0+B%*%xi%*%t(xi)%*%t(B))
      }   
      LL<-c(LL,ll)
      if(iter>5){rll<-abs(LL[length(LL)]-LL[length(LL)-1])/abs(LL[length(LL)])}
      #cat(iter,log(rll,base=10),ll," ",round(diag(S0),2)," ",round(c(B),2),"\n")
    }
    ###
  }
  list(S0=S0,B=B,ll=ll,LL=LL)
}








