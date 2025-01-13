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
num_simu = 100
n=200; p=50; q=30
lambda_seq = seq(1.0,0.1,-0.01)
alpha_seq = c(0.25,0.5,0.75)
covl_lambmin = numeric(num_simu); covl_almin = numeric(num_simu)
S_diff = numeric(num_simu); ST_diff = S_diff; lse_diff = S_diff; covl_diff = S_diff; covlpd_diff = S_diff; covreg_diff = S_diff
covl_err_est = S_diff; covl_err_pred = S_diff
lasso_TPR = numeric(num_simu); lasso_FPR = numeric(num_simu)
coef_evalmin = matrix(NA, nrow=num_simu, ncol=q+1)
Sigma_diag = matrix(NA, nrow = q+1, ncol = p); Sigma_off = matrix(NA, nrow = q+1, ncol = p*(p-1)/2)
D <- duplication.matrix(p)
Dp = solve(t(D)%*%D)%*%t(D)
diag_id = (vech(diag(p))==1)
tb_COV = matrix(NA, nrow = num_simu, ncol = 6)
Sigma_12_mat = matrix(NA, nrow=num_simu, ncol = n)
S_12_mat = Sigma_12_mat; covreg_12_mat = Sigma_12_mat; lse_12_mat = Sigma_12_mat; covl_12_mat = Sigma_12_mat

# Covariance Model Setting: 1-MA(1), 2-Clique, 3-Hub
set_model = 1
Sigma_array = array(0, dim=c(p,p,q+1))
Sigma_array[,,1] = 0.5*diag(p)
if (set_model == 1){
  Sigma_array[,,2] = sigMA(p,1,0.5)-Sigma_array[,,1]
} else if (set_model == 2){
  Sigma_array[,,2] = sigBD(p, diag0=1, a=0.5, b=0.5, k=p/5)-Sigma_array[,,1]
} else{
  Sigma_array[,,2] = sigHB(p)-Sigma_array[,,1]
}
AA = 2*Dp%*%kronecker(Sigma_array[,,1], Sigma_array[,,1])%*%t(Dp)
AB = 2*Dp%*%kronecker(Sigma_array[,,1], Sigma_array[,,2])%*%t(Dp)
BA = 2*Dp%*%kronecker(Sigma_array[,,2], Sigma_array[,,1])%*%t(Dp)
BB = 2*Dp%*%kronecker(Sigma_array[,,2], Sigma_array[,,2])%*%t(Dp)

# Covariate Setting: 1-continuous covariates, 2-binary covariates
set_covariate = 1

################# run simulations #################
for (ii in 1:num_simu){
  print(ii)
  set.seed(10000*ii)
  
  if (set_covariate == 1){
    z = matrix(runif(n*q), ncol=q)
  } else{
    z = matrix(rbinom(n*q, 1, 0.5), ncol=q)
  }
  
  Sigma_indiv = array(0, dim=c(p,p,n))#; S_indiv = Sigma_indiv; lse_indiv = Sigma_indiv; ridge_indiv = Sigma_indiv; lasso_indiv = Sigma_indiv; sgl_indiv = Sigma_indiv
  dat = matrix(0, nrow = n, ncol = p)
  for (i in 1:n){
    Sigma_indiv[,,i] = Sigma_array[,,1]+1*z[i,1]*Sigma_array[,,2] #+z[i,2]*Sigma_array[,,3]
    dat[i,] = rmvnorm(1,mean=rep(0,p),sigma=Sigma_indiv[,,i])
  }
  Sigma_coeftot = matrix(NA, nrow=q+1, ncol=p*(p+1)/2)
  for (iii in 1:(q+1)){
    Sigma_coeftot[iii,] = vech(Sigma_array[,,iii])
  }
  Sigma_coef0 = Sigma_coeftot[,diag_id]
  Sigma_coef = Sigma_coeftot[,!diag_id]
  
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
  covreg_res = covreg.em(dats,XZns,R=1,tol = 0.01)    # If we use smaller tol, the result is worse.
  for (k in 1:n){
    covreg_indiv[,,k]=covreg_res$S0+covreg_res$B%*%XZns[k,]%*%t(XZns[k,])%*%t(covreg_res$B)
  }
  print("covreg done")
  
  ####################### SparseCovReg #######################
  
  covl_cv2 = cv.spcovreg(X=XZ, Y0=Ymat_diag, Y=Ymat, lambda_seq = lambda_seq, alpha_seq = alpha_seq, k = 5, cvseed = 9999999)
  covl_lambmin[ii] = covl_cv2$lambda_min; covl_almin[ii] = covl_cv2$alpha_min
  covl_res = spcovreg(Xtilde=XZ, Ytilde0=Ymat_diag, Ytilde=Ymat, lambda1=(1-covl_almin[ii])*covl_lambmin[ii], lambda2 = covl_almin[ii]*covl_lambmin[ii], eps = 0.000001)
  print("sparsecovreg done")
  
  z_sd_mat = matrix(c(1,z_sd), nrow = q+1, ncol = p*(p-1)/2)
  z_mean_mat = matrix(c(0,z_mean), nrow = q+1, ncol = p*(p-1)/2)
  covl_beta = covl_res$beta/z_sd_mat   #  covl_res$beta/c(1,z_sd)
  covl_intadj = covl_res$beta*z_mean_mat/z_sd_mat     #  covl_res$beta*(c(0,z_mean)/c(1,z_sd))
  covl_beta[1,] = covl_beta[1,]-apply(covl_intadj,2,sum)
  z_sd_mat = matrix(c(1,z_sd), nrow = q+1, ncol = p)
  z_mean_mat = matrix(c(0,z_mean), nrow = q+1, ncol = p)
  covl_beta0 = covl_res$beta0/z_sd_mat
  covl_intadj0 = covl_res$beta0*z_mean_mat/z_sd_mat
  covl_beta0[1,] = covl_beta0[1,]-apply(covl_intadj0,2,sum)
  
  XZZ=XZ; XZZ[,-1]=z
  covl_err_est[ii] = sum(((covl_beta-Sigma_coef))^2)+sum(((covl_beta0-Sigma_coef0))^2)                  
  covl_err_pred[ii] = (sum((XZZ%*%(covl_beta-Sigma_coef))^2)+sum((XZZ%*%(covl_beta0-Sigma_coef0))^2))/n     
  
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
  
  S_indiv=array(0, c(p,p,n)); S_indiv[1:p,1:p,]=S
  ST_indiv=array(0, c(p,p,n)); ST_indiv[1:p,1:p,]=ST
  lse_indiv=array(0, c(p,p,n))
  covl_indiv=array(0, c(p,p,n))
  covlpd_indiv=array(0, c(p,p,n))
  for (i in 1:n){
    lse_indiv[,,i] = lse_coef[,,1]
    covl_indiv[,,i] = covl_coef[,,1]
    covlpd_indiv[,,i] = covlpd_coef[,,1]
    for (k in 1:q){
      lse_indiv[,,i] = lse_indiv[,,i] + z[i,k]*lse_coef[,,k+1]
      covl_indiv[,,i] = covl_indiv[,,i] + z[i,k]*covl_coef[,,k+1]
      covlpd_indiv[,,i] = covlpd_indiv[,,i] + z[i,k]*covlpd_coef[,,k+1]
    }
  }
  
  ################ Data for Figure 1 ######################
  Sigma_12_mat[ii,] = Sigma_indiv[1,2,]
  S_12_mat[ii,] = S_indiv[1,2,]
  covreg_12_mat[ii,] = covreg_indiv[1,2,]
  lse_12_mat[ii,] = lse_indiv[1,2,]
  covl_12_mat[ii,] = covl_indiv[1,2,]
  
  ####################### summary #######################
  lasso_TPR[ii] = sum(covl_coef*Sigma_array!=0)/sum(Sigma_array!=0)
  lasso_FPR[ii] = sum((covl_coef!=0)*(Sigma_array==0))/sum(Sigma_array==0)
  
  S_diff[ii] = mean(sqrt(apply((S_indiv-Sigma_indiv)^2, 3, sum)))  
  ST_diff[ii] = mean(sqrt(apply((ST_indiv-Sigma_indiv)^2, 3, sum))) 
  lse_diff[ii] = mean(sqrt(apply((lse_indiv-Sigma_indiv)^2, 3, sum)))   
  covreg_diff[ii] = mean(sqrt(apply((covreg_indiv-Sigma_indiv)^2, 3, sum))) 
  covl_diff[ii] = mean(sqrt(apply((covl_indiv-Sigma_indiv)^2, 3, sum))) 
  covlpd_diff[ii] = mean(sqrt(apply((covlpd_indiv-Sigma_indiv)^2, 3, sum)))   
  
  ####################### inference #######################
  for (k in 1:(q+1)){
    for (j in 1:p){
      Sigma_diag[k,j]=Sigma_array[j,j,k]
    }
  }
  tempind = 0
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      tempind = tempind+1
      Sigma_off[,tempind]=Sigma_array[i,j,]
    }
  }
  XZns = cbind(matrix(1,nrow=n,ncol=1), z)
  M = matrix(NA, nrow = q+1, ncol = q+1)
  for (jj in 1:(q+1)){
    xm_max = n^(1/2)+1; max_try = 10
    while (xm_max > n^(1/2)){
      M[jj,] = debiasingMatrix(t(XZns)%*%XZns/n, FALSE, n, jj, bound = 2*sqrt(log((q+1)*p*(p+1)/2)/n), linesearch = T, max_try = max_try)
      xm_max = max(abs(XZns%*%M[jj,]))
      max_try = max_try-1
    }
  }
  covl_betaU = covl_beta + (1/n)*M%*%crossprod(XZns, Ymat - XZns%*%covl_beta)
  covl_beta0U = covl_beta0 + (1/n)*M%*%crossprod(XZns, Ymat_diag - XZns%*%covl_beta0)
  covl_betaU_975 = covl_betaU*0; covl_betaU_025 = covl_betaU*0
  covl_beta0U_975 = covl_beta0U*0; covl_beta0U_025 = covl_beta0U*0
  Ahalf = M%*%t(XZns)
  Emat = Ymat - XZns%*%covl_betaU
  Emat_diag = Ymat_diag - XZns%*%covl_beta0U
  coef_err = Ahalf%*%Emat/n
  coef_err_diag = Ahalf%*%Emat_diag/n
  coef_var = covl_betaU*0
  for (iii in 1:dim(covl_betaU)[2]){
    for (jjj in 1:(q+1)){
      for (kkk in 1:n){
        coef_var[jjj,iii] = coef_var[jjj,iii] +  (Ahalf[jjj,kkk]*Emat[kkk,iii] - coef_err[jjj,iii])^2
      }
    }
  }
  coef_var = coef_var/n
  moe = 1.96/sqrt(n)*sqrt(coef_var)
  covl_betaU_975 = covl_betaU + moe
  covl_betaU_025 = covl_betaU - moe
  coef_var_diag = covl_beta0U*0
  for (iii in 1:dim(covl_beta0U)[2]){
    for (jjj in 1:(q+1)){
      for (kkk in 1:n){
        coef_var_diag[jjj,iii] = coef_var_diag[jjj,iii] +  (Ahalf[jjj,kkk]*Emat_diag[kkk,iii] - coef_err_diag[jjj,iii])^2
      }
    }
  }
  coef_var_diag = coef_var_diag/n
  moe = 1.96/sqrt(n)*sqrt(coef_var_diag)
  covl_beta0U_975 = covl_beta0U + moe
  covl_beta0U_025 = covl_beta0U - moe
  
  Vmat = matrix(NA, nrow = n, ncol = p*(p-1)/2)
  for (kk in 1:n){
    Vmat[kk,] = diag(AA+z[kk,1]*AB+z[kk,1]*BA+z[kk,1]^2*BB)[!diag_id]
  }
  MX2 = (M%*%t(XZns))^2
  moe = 1.96/n*sqrt(MX2%*%Vmat)
  true_betaU_975 = covl_betaU + moe
  true_betaU_025 = covl_betaU - moe
  
  tb_COV[ii,] = c(sum((covl_betaU_975>=Sigma_off)*(covl_betaU_025<=Sigma_off))/(prod(dim(Sigma_off))),  # total coverage
                  sum(((covl_betaU_975>=Sigma_off)*(covl_betaU_025<=Sigma_off))*(Sigma_off!=0))/sum(Sigma_off!=0),  # active coverage
                  sum(((covl_betaU_975>=Sigma_off)*(covl_betaU_025<=Sigma_off))*(Sigma_off==0))/sum(Sigma_off==0),
                  sum((true_betaU_975>=Sigma_off)*(true_betaU_025<=Sigma_off))/(prod(dim(Sigma_off))),  # total coverage
                  sum(((true_betaU_975>=Sigma_off)*(true_betaU_025<=Sigma_off))*(Sigma_off!=0))/sum(Sigma_off!=0),  # active coverage
                  sum(((true_betaU_975>=Sigma_off)*(true_betaU_025<=Sigma_off))*(Sigma_off==0))/sum(Sigma_off==0))
  
}
rm(list=c("AA","AB","BA","BB","D","Dp","covl_res","covreg_res","eigres","Emat","soft_cv","Vmat","Ymat","covl_indiv","covlpd_indiv","covreg_indiv","lse_indiv","S_indiv","Sigma_indiv","ST_indiv"))

################# save model fit result #################
# MA(1) model under Setting 1 when n=200, q=30 (adjust n and q below) 
save.image(file=paste0(path,"data/fit_n200p50q30_MA_set1_test.RData"))
# MA(1) model under Setting 2 when n=200, q=30 (adjust n and q below) 
#save.image(file=paste0(path,"data/fit_n200p50q30_MA_set2_test.RData"))
# Clique model under Setting 1 when n=200, q=30 (adjust n and q below)
#save.image(file=paste0(path,"data/fit_n200p50q30_BD_set1_test.RData"))
# Clique model under Setting 2 when n=200, q=30 (adjust n and q below)
#save.image(file=paste0(path,"data/fit_n200p50q30_BD_set2_test.RData"))
# Hub model under Setting 1 when n=200, q=30 (adjust n and q below)
#save.image(file=paste0(path,"data/fit_n200p50q30_HB_set1_test.RData"))
# Hub model under Setting 2 when n=200, q=30 (adjust n and q below)
#save.image(file=paste0(path,"data/fit_n200p50q30_HB_set2_test.RData"))
