rm(list=ls(all=TRUE))
################# set working directory #################
# please specify the path
path="/Users/rakheon_kim/Desktop/Research/SparseCovReg/"

################# load required libraries #################
library(mvtnorm)
library(matrixcalc)
library(ggplot2)
library(FinCovRegularization)
library(selectiveInference)
source(paste0(path,"code/functions/spcovreg_basic.R")) 
source(paste0(path,"code/functions/spcovreg.R"))
source(paste0(path,"code/functions/covreg_em.r")) 

# Basic Setting for n=200 and q=30 (For different n and q, adjust them below)
num_simu = 100                                   # number of simulations
n=200; p=50; q=30                                # n: sample size, p: number of response variables, q: number of covariates
lambda_seq = seq(1.0,0.1,-0.01)                  # sequence of lambda considered 
alpha_seq = c(0.25,0.5,0.75)                     # sequence of alpha considered 
covl_lambmin = numeric(num_simu); covl_almin = numeric(num_simu)   # list of selected tuning parameters, lambda and alpha, by cross validation

# Covariance Model Setting: 1-MA(1), 2-Clique, 3-Hub
set_model = 1
Sigma_array = array(0, dim=c(p,p,q+1))           # matrices for the true population-level covariance matrix and the true covariate effects
Sigma_array[,,1] = 0.5*diag(p)                   # population-level covariance matrix
if (set_model == 1){
  Sigma_array[,,2] = sigMA(p,1,0.5)-Sigma_array[,,1]       # covariate effect of the first covariate (other covariate effects are zero) 
} else if (set_model == 2){
  Sigma_array[,,2] = sigBD(p, diag0=1, a=0.5, b=0.5, k=p/5)-Sigma_array[,,1]
} else{
  Sigma_array[,,2] = sigHB(p)-Sigma_array[,,1]
}

# Covariate Setting: 1-continuous covariates, 2-binary covariates
set_covariate = 1

################# run simulations #################
for (ii in 1:num_simu){
  print(ii)
  set.seed(10000*ii)
  
  if (set_covariate == 1){
    z = matrix(runif(n*q), ncol=q)               # continous (uniform) covariates
  } else{
    z = matrix(rbinom(n*q, 1, 0.5), ncol=q)      # binary covariates
  }
  
  Sigma_indiv = array(0, dim=c(p,p,n))           # list of subject-specific covariance matrices
  dat = matrix(0, nrow = n, ncol = p)
  for (i in 1:n){
    Sigma_indiv[,,i] = Sigma_array[,,1]+1*z[i,1]*Sigma_array[,,2]
    dat[i,] = rmvnorm(1,mean=rep(0,p),sigma=Sigma_indiv[,,i])       # generate the normal response for each subject using the subject-specific covariance matrix
  }
  
  ####################### DenseSample & SparseSample #######################
  
  S=cov(dat)*(n-1)/n                                                # sample covariance matrix
  print("sample est done")
  soft_cv = threshold.cv(dat, method = "soft", thresh.len = 200, n.cv = 10, norm = "F", seed = 123)
  ST = soft(S, soft_cv$parameter.opt); diag(ST)=diag(S)             # soft thresholding estimator
  print("soft-threshold est done")
  
  ####################### DenseCovReg #######################
  
  z_mean = apply(z,2,mean); z_sd = apply(z,2,sd)
  zs = scale(z)
  XZ = cbind(matrix(1,nrow=n,ncol=1), zs)                    # scaled covariates with intercept term (colum of ones)
  Ymat = matrix(NA, nrow = n, ncol = p*(p-1)/2)              # cross-products of the demeaned response variables
  Ymat_diag = matrix(NA, nrow = n, ncol = p)                 # square of the demeaned response variables
  ij = 0
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      ij = ij + 1
      Y = (dat[,i]-mean(dat[,i]))*(dat[,j]-mean(dat[,j]))
      Ymat[,ij] = Y
    }
  }
  for (i in 1:p){
    Ymat_diag[,i]=(dat[,i]-mean(dat[,i]))^2
  }
  
  XZns = cbind(matrix(1,nrow=n,ncol=1), z)                      # non-scaled covariates with intercept term (colum of ones)
  lse_est = solve(t(XZns)%*%XZns)%*%t(XZns)%*%Ymat              # LSE estimator of the proposed covariance regression coefficients for off-diagonals
  lse_est_diag = solve(t(XZns)%*%XZns)%*%t(XZns)%*%Ymat_diag    # LSE estimator of the proposed covariance regression coefficients for diagonals
  
  ####################### CovReg #######################
  
  covreg_indiv=array(0, c(p,p,n))
  dats = scale(dat, scale=FALSE)
  covreg_res = covreg.em(dats,XZns,R=1,tol = 0.01)    # If we use smaller tol, the result is worse.
  for (k in 1:n){
    covreg_indiv[,,k]=covreg_res$S0+covreg_res$B%*%XZns[k,]%*%t(XZns[k,])%*%t(covreg_res$B)    # list of estimated subject-specific covariance matrices by Hoff and Niu (2012)
  }
  print("covreg done")
  
  ####################### SparseCovReg #######################
  
  covl_cv2 = cv.spcovreg(X=XZ, Y0=Ymat_diag, Y=Ymat, lambda_seq = lambda_seq, alpha_seq = alpha_seq, k = 5, cvseed = 9999999)
  covl_lambmin[ii] = covl_cv2$lambda_min; covl_almin[ii] = covl_cv2$alpha_min
  
}
rm(list = setdiff(ls(), c("covl_lambmin", "covl_almin")))


################# save selected tuning parameters #################
# MA(1) model under Setting 1 when n=200, q=30 (adjust n and q below) 
#save.image(file=paste0(path,"data/main/CVresult_n200p50q30_MA_set1.RData"))
# MA(1) model under Setting 2 when n=200, q=30 (adjust n and q below) 
#save.image(file=paste0(path,"data/main/CVresult_n200p50q30_MA_set2.RData"))
# Clique model under Setting 1 when n=200, q=30 (adjust n and q below)
#save.image(file=paste0(path,"data/main/CVresult_n200p50q30_BD_set1.RData"))
# Clique model under Setting 2 when n=200, q=30 (adjust n and q below)
#save.image(file=paste0(path,"data/main/CVresult_n200p50q30_BD_set2.RData"))
# Hub model under Setting 1 when n=200, q=30 (adjust n and q below)
#save.image(file=paste0(path,"data/main/CVresult_n200p50q30_HB_set1.RData"))
# Hub model under Setting 2 when n=200, q=30 (adjust n and q below)
#save.image(file=paste0(path,"data/main/CVresult_n200p50q30_HB_set2.RData"))
