rm(list=ls(all=TRUE))
################# set working directory #################
# please specify the path
path="/Users/rakheon_kim/Desktop/Research/SparseCovReg/"

################# load required libraries #################
library(mvtnorm)
library(matrixcalc)
library(ggplot2)
library(reshape2)
library(FinCovRegularization)
library(selectiveInference)
library(glmnet)
library(network)
source(paste0(path,"code/functions/spcovreg_basic.R")) 
source(paste0(path,"code/functions/spcovreg.R"))
source(paste0(path,"code/functions/covreg_em.r")) 

load(file=paste0(path,"data/main/gbmdata.RData"))
dat = genedat; z = covdat

set.seed(1)
meancv = cv.glmnet(x=z,y=dat,family = "mgaussian")
meanreg = glmnet(x=z,y=dat,family = "mgaussian", lambda = meancv$lambda.min)
dat_resid = dat - matrix(meanreg$a0, nrow = nrow(dat), ncol = ncol(dat), byrow = T)
for (kk in 1:length(meanreg$beta)){
  dat_resid[,kk] = dat_resid[,kk] - as.numeric(z%*%meanreg$beta[[kk]])
}
dat_resid=as.matrix(dat_resid)
dat = dat_resid
dat = scale(dat)

####################### estimation #######################

n = dim(dat)[1]; p = dim(dat)[2]; q = dim(z)[2]
lambda_seq = seq(1.0,0.1,-0.01)
alpha_seq = c(0.25,0.5,0.75)

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

S=cov(dat)*(n-1)/n
soft_cv = threshold.cv(dat, method = "soft", thresh.len = 200, n.cv = 10, norm = "F", seed = 123)
ST = soft(S, soft_cv$parameter.opt); diag(ST)=diag(S)
print("common est done")

covl_cv2 = cv.spcovreg(X=XZ, Y0=Ymat_diag, Y=Ymat, lambda_seq = lambda_seq, alpha_seq = alpha_seq, k = 5, cvseed = 9999999)
covl_res = spcovreg(Xtilde=XZ, Ytilde0=Ymat_diag, Ytilde=Ymat, lambda1=(1-covl_cv2$alpha_min)*covl_cv2$lambda_min, lambda2 = covl_cv2$alpha_min*covl_cv2$lambda_min, eps = 0.000001)

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

coef_evalmin = numeric(q+1)
lse_coef=array(0, c(p,p,q+1))
covl_coef=array(0, c(p,p,q+1))
covlpd_coef=array(0, c(p,p,q+1))
for (k in 1:(q+1)){
  lse_coef[,,k] = intomat(lse_est[k,], p)
  diag(lse_coef[,,k]) = lse_est_diag[k,]
  covl_coef[,,k] = intomat(covl_beta[k,], p)
  diag(covl_coef[,,k]) = covl_beta0[k,]
  coef_evalmin[k] = eigen(covl_coef[,,k])$values[p]
}
covlpd_coef=covl_coef
ccc = max(0, (-sum(coef_evalmin[-1]*(coef_evalmin[-1]<0)))-coef_evalmin[1])
epsilon = ccc/(1+ccc)
covlpd_coef = covlpd_coef*(1-epsilon)
covlpd_coef[,,1] = covlpd_coef[,,1] + epsilon*diag(p)


####################### inference #######################
XZns = cbind(matrix(1,nrow=n,ncol=1), z)
M = matrix(NA, nrow = q+1, ncol = q+1)
for (jj in 1:(q+1)){
  xm_max = n^(1/2)+1; max_try = 10
  while ((xm_max > n^(1/2))&(max_try>1)){
    M[jj,] = debiasingMatrix(t(XZns)%*%XZns/n, FALSE, n, jj, bound = 2*sqrt(log((q+1)*p*(p+1)/2)/n), linesearch = T, max_try = max_try)
    xm_max = max(abs(XZns%*%M[jj,])) # apply(abs(XZ%*%t(M)), 2, max)
    max_try = max_try-1
  }
  print(max_try)
}
covl_betaU = covl_beta + (1/n)*M%*%crossprod(XZns, Ymat - XZns%*%covl_beta)
covl_beta0U = covl_beta0 + (1/n)*M%*%crossprod(XZns, Ymat_diag - XZns%*%covl_beta0)
covl_betaU_975 = covl_betaU*0; covl_betaU_025 = covl_betaU*0
Ahalf = M%*%t(XZns)
Emat = Ymat - XZns%*%covl_betaU
coef_err = Ahalf%*%Emat/n
coef_var = covl_betaU*0
for (iii in 1:dim(covl_betaU)[2]){
  for (jjj in 1:(q+1)){
    for (kkk in 1:n){
      coef_var[jjj,iii] = coef_var[jjj,iii] +  (Ahalf[jjj,kkk]*Emat[kkk,iii] - coef_err[jjj,iii])^2
    }
  }
}
coef_var = coef_var/n
moe = 4.475862/sqrt(n)*sqrt(coef_var)        #qnorm(1-0.01/(p*(p-1)/2))
covl_betaU_975 = covl_betaU + moe
covl_betaU_025 = covl_betaU - moe
covlU_coef=array(0, c(p,p,q+1))
covlUp_coef=array(0, c(p,p,q+1))
covlLo_coef=array(0, c(p,p,q+1))
nz_coef=array(0, c(p,p,q+1))
for (k in 1:(q+1)){
  covlU_coef[,,k] = intomat(covl_betaU[k,], p)
  covlUp_coef[,,k] = intomat(covl_betaU_975[k,], p)
  covlLo_coef[,,k] = intomat(covl_betaU_025[k,], p)
  nz_coef[,,k] = (covlUp_coef[,,k]*covlLo_coef[,,k]>0)*(covl_coef[,,k]!=0)
}

rm(list=c("covl_res","Emat","Ymat","lse_coef","covlLo_coef","covlpd_coef","covlU_coef","covlUp_coef","coef_err","coef_var","covl_beta","covl_betaU","covl_betaU_025","covl_betaU_975","meancv","soft_cv"))
#save.image(file=paste0(path,"data/main/gbmdata_run.RData"))

###################### Figure 3 ######################
idsel=c( 64:66,58,61,62,57,72, 59,56,55, 25:27,15,13,17,22,16,20,21,23,24,12,14,18,19,48,49,50,45,46,47, 51:53,54,8,9,28,29,30,31,1, 40,44,43,42,41,36,37,39, 35,38,34,33,32,3,67,68,69,11,10,2,71,6,7,4,5, 73,70,63,60)
# left panel
Rn= ST[idsel,idsel]
rownames(Rn)=colnames(dat[,idsel]); colnames(Rn)=rev(colnames(dat[,idsel])); melted_R <- melt(Rn, na.rm = TRUE)
ggplot(data = melted_R, aes(Var1, rev(Var2), fill = value))+geom_tile(color = "white")+
  scale_fill_gradient2(low = "#0072B2", high = "#D55E00", mid = "white", midpoint = 0, limit = c(-0.7,0.7), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+scale_x_discrete(position = "top") +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle=90,hjust=0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank(),
        legend.position="none") + coord_fixed()
# right panel
Rn= covl_coef[idsel,idsel,1]
rownames(Rn)=colnames(dat[,idsel]); colnames(Rn)=rev(colnames(dat[,idsel])); melted_R <- melt(Rn, na.rm = TRUE)
ggplot(data = melted_R, aes(Var1, rev(Var2), fill = value))+geom_tile(color = "white")+
  scale_fill_gradient2(low = "#0072B2", high = "#D55E00", mid = "white", midpoint = 0, limit = c(-0.7,0.7), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+scale_x_discrete(position = "top") +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle=90,hjust=0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank(),
        legend.position="none") + coord_fixed()

###################### Figure 4 ######################
# left panel
Rn= covl_coef[idsel,idsel,18]+diag(p)
rownames(Rn)=colnames(dat[,idsel]); colnames(Rn)=rev(colnames(dat[,idsel])); melted_R <- melt(Rn, na.rm = TRUE)
ggplot(data = melted_R, aes(Var1, rev(Var2), fill = value))+geom_tile(color = "white")+
  scale_fill_gradient2(low = "#0072B2", high = "#D55E00", mid = "white", midpoint = 0, limit = c(-0.7,0.7), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+scale_x_discrete(position = "top") +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle=90,hjust=0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank(),
        legend.position="none") + coord_fixed()
# right panel
Rn= covl_coef[idsel,idsel,3]+diag(p)
rownames(Rn)=colnames(dat[,idsel]); colnames(Rn)=rev(colnames(dat[,idsel])); melted_R <- melt(Rn, na.rm = TRUE)
ggplot(data = melted_R, aes(Var1, rev(Var2), fill = value))+geom_tile(color = "white")+
  scale_fill_gradient2(low = "#0072B2", high = "#D55E00", mid = "white", midpoint = 0, limit = c(-0.7,0.7), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+scale_x_discrete(position = "top") +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle=90,hjust=0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank(),
        legend.position="none") + coord_fixed()

###################### Figure S1 ######################
# top left panel
Rn= covl_coef[idsel,idsel,94]+diag(p)
rownames(Rn)=colnames(dat[,idsel]); colnames(Rn)=rev(colnames(dat[,idsel])); melted_R <- melt(Rn, na.rm = TRUE)
ggplot(data = melted_R, aes(Var1, rev(Var2), fill = value))+geom_tile(color = "white")+
  scale_fill_gradient2(low = "#0072B2", high = "#D55E00", mid = "white", midpoint = 0, limit = c(-0.4,0.4), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+scale_x_discrete(position = "top") +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle=90,hjust=0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank(),
        legend.position="none") + coord_fixed()
# top right panel
Rn= covl_coef[idsel,idsel,33]+diag(p)
rownames(Rn)=colnames(dat[,idsel]); colnames(Rn)=rev(colnames(dat[,idsel])); melted_R <- melt(Rn, na.rm = TRUE)
ggplot(data = melted_R, aes(Var1, rev(Var2), fill = value))+geom_tile(color = "white")+
  scale_fill_gradient2(low = "#0072B2", high = "#D55E00", mid = "white", midpoint = 0, limit = c(-0.4,0.4), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+scale_x_discrete(position = "top") +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle=90,hjust=0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank(),
        legend.position="none") + coord_fixed()
# bottom left panel
Rn= covl_coef[idsel,idsel,28]+diag(p)
rownames(Rn)=colnames(dat[,idsel]); colnames(Rn)=rev(colnames(dat[,idsel])); melted_R <- melt(Rn, na.rm = TRUE)
ggplot(data = melted_R, aes(Var1, rev(Var2), fill = value))+geom_tile(color = "white")+
  scale_fill_gradient2(low = "#0072B2", high = "#D55E00", mid = "white", midpoint = 0, limit = c(-0.1,0.1), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+scale_x_discrete(position = "top") +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle=90,hjust=0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank(),
        legend.position="none") + coord_fixed()
# bottom right panel
Rn= covl_coef[idsel,idsel,63]+diag(p)
rownames(Rn)=colnames(dat[,idsel]); colnames(Rn)=rev(colnames(dat[,idsel])); melted_R <- melt(Rn, na.rm = TRUE)
ggplot(data = melted_R, aes(Var1, rev(Var2), fill = value))+geom_tile(color = "white")+
  scale_fill_gradient2(low = "#0072B2", high = "#D55E00", mid = "white", midpoint = 0, limit = c(-0.1,0.1), space = "Lab", name="Pearson\nCorrelation") +
  theme_minimal()+scale_x_discrete(position = "top") +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size=6),
        axis.text.x = element_text(size=6, angle=90,hjust=0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_blank(),
        legend.position="none") + coord_fixed()

###################### Figure 5 ######################
# left panel
mat = covl_coef[,,18]*nz_coef[,,18]
colnames(mat) = colnames(dat); rownames(mat) = colnames(dat)
matcol = (mat<0)*4 + (mat>0)*2
nn = as.network(mat, directed = FALSE)#, vertex.attr = mat)
set.seed(6)  # for XZ[,18]
lp = list(niter=NULL, max.delta=NULL, area=NULL, cool.exp=NULL, repulse.rad=NULL,
          init=NULL, groups=NULL, rotation=NULL, layout.control=0.5, constraints=NULL,
          round = TRUE, digits = 10)
lc = network.layout.fruchtermanreingold(nn, layout.par=lp)
par(mar = c(1, 1, 1, 1))
plot.network(nn, label = network.vertex.names(nn), coord = lc, label.cex = 0.6, vertex.cex = 0.5, edge.col = matcol, edge.lwd = abs(mat)*10)
# right panel
mat = covl_coef[,,3]*nz_coef[,,3]
colnames(mat) = colnames(dat); rownames(mat) = colnames(dat)
matcol = (mat<0)*4 + (mat>0)*2
nn = as.network(mat, directed = FALSE)#, vertex.attr = mat)
set.seed(17)  # for XZ[,3]
lp = list(niter=NULL, max.delta=NULL, area=NULL, cool.exp=NULL, repulse.rad=NULL,
          init=NULL, groups=NULL, rotation=NULL, layout.control=0.5, constraints=NULL,
          round = TRUE, digits = 10)
lc = network.layout.fruchtermanreingold(nn, layout.par=lp)
par(mar = c(1, 1, 1, 1))
plot.network(nn, label = network.vertex.names(nn), coord = lc, label.cex = 0.6, vertex.cex = 0.5, edge.col = matcol, edge.lwd = abs(mat)*10)

