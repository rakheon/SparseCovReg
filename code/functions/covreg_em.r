# The functions for implementing Hoff and Niu (2012).
# Source: http://www2.stat.duke.edu/~pdh10/Code/hoff_niu_2012_ss/covreg_em.r

# log density function for multivariate normal
# INPUT:
#  @y: p-dimensional vector of the response variable
#  @mu: p-dimensional mean vector
#  @Sig: covariance matrix
# OUTPUT: log density evaluated at y
ldmvnorm<-function(y,mu=rep(0,length(y)),Sig=diag(1,length(y)))
{
  -.5*( length(y)*log(2*pi) + log(det(Sig)) + t(y-mu)%*%solve(Sig)%*%(y-mu)  )
}

# EM algorithm function for covariance regression
# INPUT:
#  @Y: the n-by-p matrix for the p response variables
#  @X: the n-by-(q+1) design matrix for q covariates (including a column of ones for the intercept term)
#  @R: rank of the deviation between the covariance matrix and the baseline covariance matrix
#  @tol: convergence tolerance
#  @itmax: maximum number of EM iterations
# OUTPUT:
#  @S0: estimated baseline covariance matrix
#  @B: estimated regression coefficients for the latent structure
#  @ll: final log-likelihood value
#  @LL: log-likelihood trace across iterations
covreg.em<-function(Y,X,R=dim(Y)[2],tol=1e-10,itmax=1000)
{

  p<-dim(Y)[2]    # Number of response variables
  q<-dim(X)[2]    # Number of covariates 
  n<-dim(Y)[1]    # Number of observations
  S0<-cov(Y)      # Initial baseline covariance matrix (p x p)     
  iS0<-solve(S0)  # Inverse of the baseline covariance matrix
  B<-matrix(rnorm(p*q*R),p,q*R)*1e-4      # Initial estimate of regression coefficients
  iter<-0         # EM iteration counter
  LL<-NULL        # Log-likelihood trace
  rll<-10         # Relative log-likelihood change
  while( rll > tol  & iter<itmax)
  {
    B0<-B ; iter<-iter+1

    ### find expectation, var of z
    ### E-step: compute expected values of latent variables Z and their covariances
    Vz<-array(dim=c(R,R,n))       # Covariance of Z given Y and X for each subject
    Mz<-matrix(nrow=n,ncol=R)     # Mean of Z given Y and X for each subject
    for(i in 1:n)
    {
      Bx<-apply(array(B,dim=c(p,q,R)),3,"%*%",X[i,])
      Vz[,,i]<-solve(  t(Bx)%*%iS0%*%Bx + diag(R) )
      Mz[i,]<-Vz[,,i]%*%t(Bx)%*%iS0%*%Y[i,]
    }
    ###

    ### obtain MLEs
    ### M-step: Update B and baseline covariance matrix S0
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

    # Update B via least squares
    B<-t(YA)%*%XA%*%solve(t(XA)%*%XA) 
    # Update baseline covariance
    E<-YA-XA%*%t(B)
    S0<- (t(E)%*%E)/dim(Y)[1]  ; iS0<-solve(S0)
    ###


    ### Compute log-likelihood
    if(iter%%5==0)
    {
      ll<-0
      for(i in 1:dim(Y)[1])
      {
        xi<-matrix(outer(X[i,],diag(R)),nrow=R*q,ncol=R)
        ll<-ll+ldmvnorm(Y[i,],Sig=S0+B%*%xi%*%t(xi)%*%t(B))
      }   
      LL<-c(LL,ll)
      if(iter>5){rll<-abs(LL[length(LL)]-LL[length(LL)-1])/abs(LL[length(LL)])}  # Relative log-likelihood change 
      #cat(iter,log(rll,base=10),ll," ",round(diag(S0),2)," ",round(c(B),2),"\n")
    }
    ###
  }
  list(S0=S0,B=B,ll=ll,LL=LL)
}








