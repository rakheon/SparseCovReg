# objective function
# INPUT:
#  @Xtilde: the n-by-(q+1) design matrix for q covariates (including a column of ones for the intercept term)
#  @Ytilde0: the n-by-p matrix for the p squared zero-mean response variables
#  @Ytilde: the n-by-(p*(p-1)/2) matrix for the cross-product of p zero-mean response variables
#  @beta0: the (q+1)-by-p matrix of the estimator of diagonal entries in B_0,B_1,...,B_q 
#  @beta: the (q+1)-by-(p*(p-1)/2) matrix of the estimator of off-diagonal entries in B_0,B_1,...,B_q
#  @lambda1: tuning parameter for the group lasso penalty
#  @lambda2: tuning parameter for the lasso penalty
# OUTPUT: the objective function value evaluated at beta0 and beta
covlasso <- function(Xtilde, Ytilde0, Ytilde, beta0, beta, lambda1, lambda2){
  return(1/(2*dim(Xtilde)[1])*(crossprod(c(Ytilde0 - Xtilde %*% beta0))+crossprod(c(Ytilde - Xtilde %*% beta))) + lambda1*(sum(sqrt(rowSums(cbind(beta0[-1,], beta[-1,])^2)))) + lambda2*(sum(abs(beta)) + sum(abs(beta0[-1,]))))
}

# spcovreg function
# INPUT:
#  @Xtilde: the n-by-(q+1) design matrix for q covariates (including a column of ones for the intercept term)
#  @Ytilde0: the n-by-p matrix for the p squared zero-mean response variables
#  @Ytilde: the n-by-(p*(p-1)/2) matrix for the cross-product of p zero-mean response variables
#  @lambda1: tuning parameter for the group lasso penalty
#  @lambda2: tuning parameter for the lasso penalty
#  @beta_start: Can be used if a specific initial estimator of beta is preferred
#  @beta0_start: Can be used if a specific initial estimator of beta0 is preferred
#  @eps: convergence tolerance
# OUTPUT:
#  @beta0: the (q+1)-by-p matrix of the estimator of diagonal entries in B_0,B_1,...,B_q 
#  @beta: the (q+1)-by-(p*(p-1)/2) matrix of the estimator of off-diagonal entries in B_0,B_1,...,B_q
#  @fmin: a vector of objective function values for all iterations
#  @iter: number of iterations
spcovreg <- function(Xtilde, Ytilde0, Ytilde, lambda1, lambda2, beta_start = NULL, beta0_start = NULL, eps = 0.0001){

  # set dimension 
  n <- dim(Xtilde)[1]
  p <- dim(Xtilde)[2]
  q <- dim(Ytilde)[2]
  d <- dim(Ytilde0)[2]
  
  if(is.null(beta_start)){
    beta_start <- matrix(1, nrow=p, ncol=q)
    beta0_start <- matrix(1, nrow=p, ncol=d)
  } #else if (length(beta_start) != p){
    #stop("Number of parameters should be compatible with the number of predictors")
  #}
  # Coordinate-descent implementation. Stop when the difference between objective functions is less than eps.
  beta <- beta_start
  beta0 <- beta0_start
  # We will use FULL-residual approach
  full_res0 <- Ytilde0 - (Xtilde %*% beta0)
  full_res <- Ytilde - (Xtilde %*% beta)
  iter <- 0; itermax <- 1000
  diff <- 1
  fmin_old <- covlasso(Xtilde, Ytilde0, Ytilde, beta0_start, beta_start, lambda1, lambda2)
  fmin_vec = NULL
  scaleconst = numeric(p)
  while(diff >= eps & iter <= itermax){
    iter <- iter + 1
    beta0_old = beta0
    beta_old = beta
    beta0_temp = beta0
    beta_temp = beta
    # intercept: only L1 penalty
    beta0_temp[1,] = t(Xtilde[,1])%*%(full_res0+as.matrix(Xtilde[,1],ncol=1)%*%beta0[1,])/n
    beta_temp[1,] = t(Xtilde[,1])%*%(full_res+as.matrix(Xtilde[,1],ncol=1)%*%beta[1,])/n
    if (iter==1) scaleconst[1] = as.numeric(crossprod(Xtilde[,1])/n)
    normconst = as.numeric(1/(scaleconst[1] ))
    beta[1,] <- normconst*soft(scaleconst[1]*beta[1,]+t(Xtilde[,1])%*%full_res/n, lambda2 )  # 1/(sum(Xtilde[,i]^2)/n + lambda1/sqrt(sum(beta0[i,]^2)+sum(beta[i,]^2)) )*soft(t(Xtilde[,i])%*%partresid/n, lambda2 ) 
    beta0[1,] <- normconst*(scaleconst[1]*beta0[1,]+t(Xtilde[,1])%*%full_res0/n)  # 1/(sum(Xtilde[,i]^2)/n + lambda1/sqrt(sum(beta0[i,]^2)+sum(beta[i,]^2)) )*(t(Xtilde[,i])%*%partresid0/n)
    full_res <- full_res + as.matrix(Xtilde[,1],ncol=1)%*%(beta_old[1,]-beta[1,])
    full_res0 <- full_res0 + as.matrix(Xtilde[,1],ncol=1)%*%(beta0_old[1,]-beta0[1,])
    # slope for each covariate: L2 + L1 penalty
    for (i in 2:p){
      beta0_temp[i,] = t(Xtilde[,i])%*%(full_res0+as.matrix(Xtilde[,i],ncol=1)%*%beta0[i,])/n
      beta_temp[i,] = t(Xtilde[,i])%*%(full_res+as.matrix(Xtilde[,i],ncol=1)%*%beta[i,])/n
      if (iter==1) scaleconst[i] = as.numeric(crossprod(Xtilde[,i])/n)
      cond = (sqrt(crossprod(soft(beta0_temp[i,], lambda2) )+crossprod(soft(beta_temp[i,], lambda2) )) < lambda1) # (sqrt(sum(beta0_temp[i,]^2)+sum(beta_temp[i,]^2)) < lambda1)
      if (cond==TRUE){
        beta[i,] = beta[i,]*0
        beta0[i,] = beta0[i,]*0
        full_res <- full_res + as.matrix(Xtilde[,i],ncol=1)%*%(beta_old[i,]-beta[i,])
        full_res0 <- full_res0 + as.matrix(Xtilde[,i],ncol=1)%*%(beta0_old[i,]-beta0[i,])
      } else{
        normconst = as.numeric(1/(scaleconst[i] + lambda1/sqrt(crossprod(beta0[i,])+crossprod(beta[i,])) ))
        beta[i,] <- normconst*soft(scaleconst[i]*beta[i,]+t(Xtilde[,i])%*%full_res/n, lambda2 )  # 1/(sum(Xtilde[,i]^2)/n + lambda1/sqrt(sum(beta0[i,]^2)+sum(beta[i,]^2)) )*soft(t(Xtilde[,i])%*%partresid/n, lambda2 ) 
        beta0[i,] <- normconst*soft(scaleconst[i]*beta0[i,]+t(Xtilde[,i])%*%full_res0/n, lambda2 )  # 1/(sum(Xtilde[,i]^2)/n + lambda1/sqrt(sum(beta0[i,]^2)+sum(beta[i,]^2)) )*(t(Xtilde[,i])%*%partresid0/n)
        full_res <- full_res + as.matrix(Xtilde[,i],ncol=1)%*%(beta_old[i,]-beta[i,])
        full_res0 <- full_res0 + as.matrix(Xtilde[,i],ncol=1)%*%(beta0_old[i,]-beta0[i,])
      }
    }
    fmin <- covlasso(Xtilde, Ytilde0, Ytilde, beta0, beta, lambda1, lambda2)
    fmin_vec <- c(fmin_vec , fmin)
    diff <- fmin_old - fmin # mean(abs(beta0_temp - beta0)) # fmin_old - fmin
    #print(fmin_old - fmin)
    fmin_old <- fmin
  }
  # Return beta - the solution, fmin - optimal function value (value of objective at beta)
  return(list(beta0 = beta0, beta = beta, fmin = fmin_vec, iter = iter))
}

# cross-validation for spcovreg
# INPUT:
#  @X: the n-by-(q+1) design matrix for q covariates (including a column of ones for the intercept term)
#  @Y0: the n-by-p matrix for the p squared zero-mean response variables
#  @Y: the n-by-(p*(p-1)/2) matrix for the cross-product of p zero-mean response variables
#  @lambda_seq: sequence of values for penalty tuning parameter (lambda)
#  @alpha_seq: sequence of values for weight tuning parameter (alpha)
#  @k: number of folds in K-fold cross-validation
#  @cvseed: If fixing the seed is preferred
# OUTPUT:
#  @lambda_min: the lambda value which minimizes the cross-validation error
#  @alpha_min: the alpha value which minimizes the cross-validation error
#  @cvm: the matrix of cross-validation errors computed on a grid of of alpha and lambda
#  @cvse: the standard error of cross-validation errors computed on a grid of of alpha and lambda
cv.spcovreg <- function(X, Y0, Y, lambda_seq, alpha_seq, k = 5, cvseed=1234){
  
  set.seed(cvseed)

  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  d <- dim(Y0)[2]
  
  idfold <- sample(1:n) %% k + 1
  
  n_lambda <- length(lambda_seq)
  n_alpha <- length(alpha_seq)
  sqerror <- array(NA, c(n, n_lambda, n_alpha))
  cvfold <- array(NA, c(k, n_lambda, n_alpha))
  for (fold in 1:k){
    print(c("fold: ", fold))
    # Training data
    xtrain = X[idfold != fold, ]
    ytrain = Y[idfold != fold, ]
    y0train = Y0[idfold != fold, ]
    
    # Test data
    xtest = X[idfold == fold, ]
    ytest = Y[idfold == fold, ]
    y0test = Y0[idfold == fold, ]
    
    for (l in 1:n_lambda){
      print(l)
      for (m in 1:n_alpha){
        cvfit <- spcovreg(Xtilde=xtrain, Ytilde0=y0train, Ytilde=as.matrix(ytrain), lambda1=(1-alpha_seq[m])*lambda_seq[l], lambda2 = alpha_seq[m]*lambda_seq[l], eps = 0.000001)
        testfitted <- xtest %*% cvfit$beta
        testfitted0 <- xtest %*% cvfit$beta0
        cvfold[fold, l, m] <- mean(c(colMeans((ytest - testfitted)^2), colMeans((y0test - testfitted0)^2)))
        sqerror[idfold == fold, l, m] <- rowSums((ytest - testfitted)^2)+rowSums((y0test - testfitted0)^2)
      }
    }
    
  }
  # Calculate CV(lambda) and SE_CV(lambda) for each value of lambda
  cvm <- apply(cvfold,c(2,3), mean) # apply(sqerror,c(2,3), mean) 
  cvse <- apply(cvfold,c(2,3), sd)/sqrt(k)
  
  # Find lambda_min
  lambda_tb = matrix(NA, nrow=n_lambda*n_alpha, ncol=3)
  lambda_tb[,1] = 1:(n_lambda*n_alpha)
  lambda_tb[,2] = rep(1:n_lambda, n_alpha)
  lambda_tb[,3] = rep(1:n_alpha, each=n_lambda)
  lambda_min <- lambda_seq[lambda_tb[which.min(cvm),2]]
  alpha_min <- alpha_seq[lambda_tb[which.min(cvm),3]]
  
  return(list(lambda_min = lambda_min, alpha_min = alpha_min, cvm = cvm, cvse = cvse)) #, cvfold = cvfold, sqerror = sqerror))
}

