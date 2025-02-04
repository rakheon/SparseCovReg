# soft-thresholding
# INPUT:
#  @a: the input scalar/vector/matrix
#  @lambda: the threshold
# OUTPUT: the soft-thresheld scalar/vector/matrix
soft <- function(a, lambda){
  return(sign(a)*pmax(abs(a)-lambda, 0))
}

# inverse function of vectorization
# INPUT:
#  @vec: the vector of off-diagonal elements in a p by p matrix
#  @p: the dimension of the square matrix
# OUTPUT: the p by p matrix where off-diagonals are equal to the elements of vec
intomat = function(vec, p){
  A = matrix(0, p, p)
  ij = 0
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      ij = ij + 1
      A[i,j] = vec[ij]
    }
  }
  A = A + t(A)
  return(A)
}

# true covariance matrix generator for the block (Clique) model
# INPUT:
#  @p: the dimension of the covariance matrix
#  @diag0: the diagonals of the covariance matrix
#  @a: the lower limit of the non-zero off-diagonal elements
#  @b: the upper limit of the non-zero off-diagonal elements
#  @k: block size
#  @hetero: if the heteroskedastic diagonals are needed
# OUTPUT: the generated true covariance matrix
sigBD = function(p, diag0, a, b, k, hetero=NULL){
    # a and b are the limits of the uniform distribution and k is the block size
    A = diag(rep(diag0, p))
    m = p / k
    for (h in 1 : m){
        for (i in ((h - 1) * k + 1) : (h * k)){
            for (j in ((h - 1) * k + 1) : (h * k)){
                if (i > j){
                    temp = runif(1, a, b)
                    A[i, j] = temp
                    A[j, i] = temp
                }
            }
        }
    }
    if (!is.null(hetero)){
        D = runif(p, 0.1, 10)
        #D = runif(p, 2, 5000)
        A = diag(D) %*% A %*% diag(D)
    }
    return(A)
}

# true covariance matrix generator for the MA model
# INPUT:
#  @p: the dimension of the covariance matrix
#  @diag0: the diagonals of the covariance matrix
#  @offdiag: a vector of values where each value represents the value of off-diagonal elements 
#  @hetero: if the heteroskedastic diagonals are needed
# OUTPUT: the generated true covariance matrix
sigMA = function(p, diag0=1, offdiag, hetero=NULL){
    # a and b are the limits of the uniform distribution and k is the block size
    A = diag(rep(diag0, p))
    for (i in 1:p){
        for (j in 1:p){
            for (k in 1:length(offdiag)){
                if (abs(i-j)==k){A[i,j] <- offdiag[k]}
            }
            #if (abs(i-j)==1){A[i,j] <- a}
        }
    }
    if (!is.null(hetero)){
        D = runif(p, 0.1, 10)
        A = diag(D) %*% A %*% diag(D)
    }
    return(A)
}

# true covariance matrix generator for the Hub model
# INPUT:
#  @p: the dimension of the covariance matrix
# OUTPUT: the generated true covariance matrix as described in Kim and Zhang (2025?)
sigHB = function(p){
  Sigma <- matrix(0, nrow = p, ncol = p)
  Sigma[1, 1:(p / 10)] <- sign(rnorm(n = p / 10))
  Sigma[(p / 10 + 1), (p / 10 + 1):(2 * p / 10)] <- sign(rnorm(n = p / 10))
  Sigma[(2 * p / 10 + 1), (2 * p / 10 + 1):(3 * p / 10)] <- sign(rnorm(n = p / 10))
  Sigma[(3 * p / 10 + 1), (3 * p / 10 + 1):(4 * p / 10)] <- sign(rnorm(n = p / 10))
  Sigma[(4 * p / 10 + 1), (4 * p / 10 + 1):(5 * p / 10)] <- sign(rnorm(n = p / 10))
  Sigma[(5 * p / 10 + 1), (5 * p / 10 + 1):(6 * p / 10)] <- sign(rnorm(n = p / 10))
  Sigma[(6 * p / 10 + 1), (6 * p / 10 + 1):(7 * p / 10)] <- sign(rnorm(n = p / 10))
  Sigma[(7 * p / 10 + 1), (7 * p / 10 + 1):(8 * p / 10)] <- sign(rnorm(n = p / 10))
  Sigma[(8 * p / 10 + 1), (8 * p / 10 + 1):(9 * p / 10)] <- sign(rnorm(n = p / 10))
  Sigma[(9 * p / 10 + 1), (9 * p / 10 + 1):p] <- sign(rnorm(n = p / 10))
  Sigma=abs(Sigma)
  diag(Sigma) <- 0
  ind <- lower.tri(Sigma, diag = FALSE)
  Sigma[ind] <- t(Sigma)[ind] # symmetrize
  Sigma=0.4*Sigma; diag(Sigma)=1
  A=Sigma
  return(A)
}

