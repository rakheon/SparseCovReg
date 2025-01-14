rm(list=ls(all=TRUE))
################# set working directory #################
# please specify the path
path="/Users/rakheon_kim/Desktop/Research/SparseCovReg/"

################# load required libraries #################
library(mvtnorm)
library(matrixcalc)
library(ggplot2)
library(FinCovRegularization)
library(covreg)
source(paste0(path,"code/functions/spcovreg_basic.R")) 
source(paste0(path,"code/functions/spcovreg.R"))
source(paste0(path,"code/functions/covreg_em.r")) 

# Basic Setting for q=1 (For different q, adjust below)
num_simu = 20
n=200; p=2; q=1    # q: either 1,3,5 or 10
lambda_seq = seq(1.0,0.1,-0.01)
alpha_seq = c(0.25,0.5,0.75)
covl_lambmin = numeric(num_simu); covl_almin = numeric(num_simu)
S_diff = numeric(num_simu); ST_diff = S_diff; lse_diff = S_diff; covl_diff = S_diff; covlpd_diff = S_diff; covreg_diff = S_diff; covreg_mcmc_diff = S_diff; covl_diff3 = S_diff
coef_evalmin = matrix(NA, nrow=num_simu, ncol=q+1)

################# run simulations #################
for (ii in 1:num_simu){
  print(ii)
  set.seed(10000*ii)
  
  z = matrix(runif(n*q), ncol=q)
  if (q==1){
    z = cbind(z, z^2)
  } else if(q==3){
    z = cbind(z, z[,1]*z, z[,2]*z[,2:3], z[,3]*z[,3])
  } else if(q==5){
    z = cbind(z, z[,1]*z, z[,2]*z[,2:5], z[,3]*z[,3:5], z[,4]*z[,4:5], z[,5]*z[,5])
  } else {
    z = cbind(z, z[,1]*z, z[,2]*z[,2:10], z[,3]*z[,3:10], z[,4]*z[,4:10], z[,5]*z[,5:10], z[,6]*z[,6:10], z[,7]*z[,7:10], z[,8]*z[,8:10], z[,9]*z[,9:10], z[,10]*z[,10])
  }
  B0 = matrix(c(1,-1,1,1), ncol=2)
  C0 = B0 %*% diag(c(1,1/3)) %*% t(B0) # diag(2) # B0 %*% diag(c(1,1/3)) %*% t(B0)
  z1 = cbind(matrix(1,nrow=n,ncol=1), z[,1])
  Sigma_indiv = array(0, dim=c(p,p,n))
  dat = matrix(0, nrow = n, ncol = p)
  for (i in 1:n){
    Sigma_indiv[1:2,1:2,i] = (C0/2) + (B0/2)%*%z1[i,]%*%t(z1[i,])%*%t(B0/2)
    dat[i,] = rmvnorm(1,mean=rep(0,p),sigma=Sigma_indiv[,,i])
  }
  
  ####################### DenseSample & SparseSample #######################
  
  S=cov(dat)*(n-1)/n  
  print("sample est done")
  soft_cv = threshold.cv(dat, method = "soft", thresh.len = 200, n.cv = 10, norm = "F", seed = 123)
  ST = soft(S, soft_cv$parameter.opt); diag(ST)=diag(S)
  print("soft-threshold est done")
  
  ####################### DenseCovReg #######################
  
  z_mean = apply(z,2,mean); z_sd = apply(z,2,sd)
  zs = scale(z)
  XZ = cbind(matrix(1,nrow=n,ncol=1), zs)
  Ymat = matrix(NA, nrow = n, ncol = p*(p-1)/2)
  Ymat_diag = matrix(NA, nrow = n, ncol = p)
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
  
  XZns = cbind(matrix(1,nrow=n,ncol=1), z)
  lse_est = solve(t(XZns)%*%XZns)%*%t(XZns)%*%Ymat
  lse_est_diag = solve(t(XZns)%*%XZns)%*%t(XZns)%*%Ymat_diag
  
  ####################### CovReg #######################

  covreg_indiv=array(0, c(p,p,n))
  dats = scale(dat, scale=FALSE)
  covreg_res = covreg.em(dats,XZns[,1:(q+1)],R=1,tol = 0.000000001)
  for (k in 1:n){
    covreg_indiv[,,k]=covreg_res$S0+covreg_res$B%*%XZns[k,1:(q+1)]%*%t(XZns[k,1:(q+1)])%*%t(covreg_res$B)
  }
  covreg_indiv_mcmc=array(0, c(p,p,n))
  datfr = as.data.frame(cbind(scale(dat, scale = F, center = T), z[,1:q])) # XZ[,2:(q+1)])) # z[,1:q]))
  r_names = paste0("r", 1:p)
  c_names = paste0("c", 1:q)
  colnames(datfr)=c(r_names,c_names)
  fmean=as.formula(cbind(r1,r2)~1) #as.formula(cbind(r1,r2)~c1+c2)
  if (q==1){
    fcov=as.formula(cbind(r1,r2)~c1)
  } else if(q==3){
    fcov=as.formula(cbind(r1,r2)~c1+c2+c3)
  } else if(q==5){
    fcov=as.formula(cbind(r1,r2)~c1+c2+c3+c4+c5)
  } else {
    fcov=as.formula(cbind(r1,r2)~c1+c2+c3+c4+c5+c6+c7+c8+c9+c10)
  }
  fit=covreg.mcmc(fmean, fcov, data=datfr, R=1) #, niter = 50000, nthin = 50)
  M.psamp=m.psamp(fit)
  S.psamp=cov.psamp(fit)
  S.psamp=S.psamp[,,,-(1:200)] # burning period
  ps_res = apply(S.psamp,c(1,2,3),mean)
  for (k in 1:n){
    covreg_indiv_mcmc[,,k]=ps_res[k,,]
  }
  print("covreg done")
  
  ####################### SparseCovReg #######################
  
  covl_cv2 = cv.spcovreg(X=XZ, Y0=Ymat_diag, Y=Ymat, lambda_seq = lambda_seq, alpha_seq = alpha_seq, k = 5, cvseed = 9999999)
  covl_lambmin[ii] = covl_cv2$lambda_min; covl_almin[ii] = covl_cv2$alpha_min
  covl_res = spcovreg(Xtilde=XZ, Ytilde0=Ymat_diag, Ytilde=Ymat, lambda1=(1-covl_almin[ii])*covl_lambmin[ii], lambda2 = covl_almin[ii]*covl_lambmin[ii], eps = 0.000001)
  print("sparsecovreg done")
  
  z_sd_mat = matrix(c(1,z_sd), nrow = ncol(XZ), ncol = p*(p-1)/2)
  z_mean_mat = matrix(c(0,z_mean), nrow = ncol(XZ), ncol = p*(p-1)/2)
  covl_beta = covl_res$beta/z_sd_mat   #  covl_res$beta/c(1,z_sd)
  covl_intadj = covl_res$beta*z_mean_mat/z_sd_mat     #  covl_res$beta*(c(0,z_mean)/c(1,z_sd))
  covl_beta[1,] = covl_beta[1,]-apply(covl_intadj,2,sum)
  z_sd_mat = matrix(c(1,z_sd), nrow = ncol(XZ), ncol = p)
  z_mean_mat = matrix(c(0,z_mean), nrow = ncol(XZ), ncol = p)
  covl_beta0 = covl_res$beta0/z_sd_mat
  covl_intadj0 = covl_res$beta0*z_mean_mat/z_sd_mat
  covl_beta0[1,] = covl_beta0[1,]-apply(covl_intadj0,2,sum)
  
  ####################### SparseCovReg (linear) #######################
  
  covl_cv3 = cv.spcovreg(X=XZ[,1:(q+1)], Y0=Ymat_diag, Y=Ymat, lambda_seq = lambda_seq, alpha_seq = alpha_seq, k = 5, cvseed = 9999999)
  covl_res3 = spcovreg(Xtilde=XZ[,1:(q+1)], Ytilde0=Ymat_diag, Ytilde=Ymat, lambda1=(1-covl_cv3$alpha_min)*covl_cv3$lambda_min, lambda2 = covl_cv3$alpha_min*covl_cv3$lambda_min, eps = 0.000001)
  print("covl done")
  
  z_sd_mat = matrix(c(1,z_sd[1:q]), nrow = ncol(XZ[,1:(q+1)]), ncol = p*(p-1)/2)
  z_mean_mat = matrix(c(0,z_mean[1:q]), nrow = ncol(XZ[,1:(q+1)]), ncol = p*(p-1)/2)
  covl_beta3 = covl_res3$beta/z_sd_mat   #  covl_res$beta/c(1,z_sd)
  covl_intadj = covl_res3$beta*z_mean_mat/z_sd_mat     #  covl_res$beta*(c(0,z_mean)/c(1,z_sd))
  covl_beta3[1,] = covl_beta3[1,]-apply(covl_intadj,2,sum)
  z_sd_mat = matrix(c(1,z_sd[1:q]), nrow = ncol(XZ[,1:(q+1)]), ncol = p)
  z_mean_mat = matrix(c(0,z_mean[1:q]), nrow = ncol(XZ[,1:(q+1)]), ncol = p)
  covl_beta03 = covl_res3$beta0/z_sd_mat
  covl_intadj0 = covl_res3$beta0*z_mean_mat/z_sd_mat
  covl_beta03[1,] = covl_beta03[1,]-apply(covl_intadj0,2,sum)
  
  ####################### Saving results #######################
  
  lse_coef=array(0, c(p,p,q+1))
  covl_coef=array(0, c(p,p,q+1))
  covlpd_coef=array(0, c(p,p,q+1))
  for (k in 1:(q+1)){
    lse_coef[,,k] = intomat(lse_est[k,], p)
    diag(lse_coef[,,k]) = lse_est_diag[k,]
    covl_coef[,,k] = intomat(covl_beta[k,], p)
    diag(covl_coef[,,k]) = covl_beta0[k,]
    coef_evalmin[ii,k] = eigen(covl_coef[,,k])$values[p]
  }
  covlpd_coef=covl_coef
  tempmat = covl_coef[,,1]
  for (k in 2:(q+1)){
    if (coef_evalmin[ii,k]<0){
      eigres = eigen(covl_coef[,,k])
      coefneg = eigres$vectors%*%diag(pmin(eigres$values,0))%*%t(eigres$vectors)
      tempmat = tempmat + coefneg
    }
  }
  ccc = max(0, -eigen(tempmat)$values[p])
  epsilon = ccc/(1+ccc)
  covlpd_coef = covlpd_coef*(1-epsilon)
  covlpd_coef[,,1] = covlpd_coef[,,1] + epsilon*diag(p)
  
  covl_coef3=array(0, c(p,p,ncol(XZ[,1:(q+1)])))
  for (k in 1:ncol(XZ[,1:(q+1)])){
    covl_coef3[,,k] = intomat(covl_beta3[k,], p)
    diag(covl_coef3[,,k]) = covl_beta03[k,]
  }
  
  S_indiv=array(0, c(p,p,n)); S_indiv[1:p,1:p,]=S
  ST_indiv=array(0, c(p,p,n)); ST_indiv[1:p,1:p,]=ST
  lse_indiv=array(0, c(p,p,n))
  covl_indiv=array(0, c(p,p,n))
  covlpd_indiv=array(0, c(p,p,n))
  covl_indiv3=array(0, c(p,p,n))
  for (i in 1:n){
    lse_indiv[,,i] = lse_coef[,,1]
    covl_indiv[,,i] = covl_coef[,,1]
    covlpd_indiv[,,i] = covlpd_coef[,,1]
    covl_indiv3[,,i] = covl_coef3[,,1]
    for (k in 1:q){
      lse_indiv[,,i] = lse_indiv[,,i] + z[i,k]*lse_coef[,,k+1]
      covl_indiv[,,i] = covl_indiv[,,i] + z[i,k]*covl_coef[,,k+1]
      covlpd_indiv[,,i] = covlpd_indiv[,,i] + z[i,k]*covlpd_coef[,,k+1]
      covl_indiv3[,,i] = covl_indiv3[,,i] + z[i,k]*covl_coef3[,,k+1]
    }
  }
  
  ####################### summary #######################
  
  S_diff[ii] = mean(sqrt(apply((S_indiv-Sigma_indiv)^2, 3, sum)))  
  ST_diff[ii] = mean(sqrt(apply((ST_indiv-Sigma_indiv)^2, 3, sum))) 
  lse_diff[ii] = mean(sqrt(apply((lse_indiv-Sigma_indiv)^2, 3, sum)))   
  covreg_diff[ii] = mean(sqrt(apply((covreg_indiv-Sigma_indiv)^2, 3, sum)))
  covreg_mcmc_diff[ii] = mean(sqrt(apply((covreg_indiv_mcmc-Sigma_indiv)^2, 3, sum)))
  covl_diff[ii] = mean(sqrt(apply((covl_indiv-Sigma_indiv)^2, 3, sum))) 
  covlpd_diff[ii] = mean(sqrt(apply((covlpd_indiv-Sigma_indiv)^2, 3, sum)))   
  covl_diff3[ii] = mean(sqrt(apply((covl_indiv3-Sigma_indiv)^2, 3, sum))) 
  
}

########## Save results when q=1 (adjust the file name below for different q) ############
#save.image(file=paste0(path,"data/supp/supp_q1.RData"))


########## Numbers for Table S1 ############
mean(covreg_diff); mean(covreg_mcmc_diff); mean(covl_diff3); mean(covl_diff)
sd(covreg_diff); sd(covreg_mcmc_diff); sd(covl_diff3); sd(covl_diff)

## code for plots ##
n=200; p=2; q=3; ii=1
set.seed(10000*ii)
z = matrix(runif(n*q), ncol=q)
z = cbind(z, z[,1]*z, z[,2]*z[,2:3], z[,3]*z[,3])
B0 = matrix(c(1,-1,1,1), ncol=2)
C0 = B0 %*% diag(c(1,1/3)) %*% t(B0) # diag(2) # B0 %*% diag(c(1,1/3)) %*% t(B0)
z1 = cbind(matrix(1,nrow=n,ncol=1), z[,1])
Sigma_indiv = array(0, dim=c(p,p,n))
dat = matrix(0, nrow = n, ncol = p)
for (i in 1:n){
  Sigma_indiv[1:2,1:2,i] = (C0/2) + (B0/2)%*%z1[i,]%*%t(z1[i,])%*%t(B0/2)
  dat[i,] = rmvnorm(1,mean=rep(0,p),sigma=Sigma_indiv[,,i])
}
plot(z[,1], Sigma_indiv[1,2,], col="white", pch=16, ylim=c(-6,2), xlab=NA, ylab=NA)
points(z[,1], (dat[,1]-mean(dat[,1]))*(dat[,2]-mean(dat[,2])), cex=0.5, pch=16)
plot(z[,2], Sigma_indiv[1,2,], col="white", pch=16, ylim=c(-6,2), xlab=NA, ylab=NA)
points(z[,2], (dat[,1]-mean(dat[,1]))*(dat[,2]-mean(dat[,2])), cex=0.5, pch=16)




