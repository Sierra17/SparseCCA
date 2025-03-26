# Alternating CCA Scheme

# CCA solver for matrices X, Y with n rows and p, q columns by direct approach
# requires inverting Sxx and Syy
cca_solver <- function(X, Y, k, n){
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  
  # May be better to use t(X) %*% Y for sparse matrices
  Sxx <- crossprod(Xc) / (n - 1)
  Syy <- crossprod(Yc) / (n - 1)
  Sxy <- crossprod(Xc, Yc) / (n - 1)
  
  specSxx <- eigen(Sxx)
  invSxx_sqrt <- eigSxx$vectors %*% diag(1 / sqrt(specSxx$values)) %*% t(specSxx$vectors)
  specSyy <- eigen(Syy)
  invSyy_sqrt <- eigSyy$vectors %*% diag(1 / sqrt(specSyy$values)) %*% t(specSyy$vectors)
  
  # Form the matrix for SVD: invSxx_sqrt * Sxy * invSyy_sqrt
  M <- invSxx_sqrt %*% Sxy %*% invSyy_sqrt
  
  svdM <- svd(M)
  
  can_corrs <- svdM$d[1:k]
  
  A <- invSxx_sqrt %*% svdM$u[, 1:k]
  B <- invSyy_sqrt %*% svdM$v[, 1:k]
  
  return(list(canonical_correlations = can_corrs, x_weights = A, y_weights = B))
}

# Alternating approach to CCA for matrices X, Y with n rows and p, q columns
# matrices A, B have to be initialized separately

alternating_cca_solver <- function(X, Y, k, n, A, B, tol = 1e-6, max.iter = 100){
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  
  Sxx <- crossprod(Xc) / (n - 1)
  Syy <- crossprod(Yc) / (n - 1)
  
  for (iter in 1:max_iter) {
    # --- Update A given B ---
    # Solve least squares: A = (X'X)^{-1} X'Y B
    A <- solve(Sxx, crossprod(Xc, Yc) %*% B)
    # Normalize A so that A' Sxx A = I
    norm_A <- sqrt(diag(t(A) %*% Sxx %*% A))
    A <- A %*% diag(1 / norm_A)
    
    # --- Update B given A ---
    # Solve least squares: B = (Y'Y)^{-1} Y'X A
    B <- solve(Syy, crossprod(Yc, Xc) %*% A)
    # Normalize B so that B' Syy B = I
    norm_B <- sqrt(diag(t(B) %*% Syy %*% B))
    B <- B %*% diag(1 / norm_B)
    
    # Check convergence (using the maximum change in A)
    if (max(abs(A - A_old)) < tol) break
    A_old <- A
  }
}



cca_alternating <- function(X, Y, k, tol = 1e-6, max_iter = 100) {
  # Center the data
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  n <- nrow(Xc)
  
  # Dimensions of X and Y
  p <- ncol(Xc)
  q <- ncol(Yc)
  
  # Sample covariance matrices (using denominator n-1)
  Sxx <- crossprod(Xc) / (n - 1)  # p x p
  Syy <- crossprod(Yc) / (n - 1)  # q x q
  
  # Initial B: a random q x k matrix with orthonormal columns
  set.seed(123)  # For reproducibility (optional)
  B <- matrix(rnorm(q * k), nrow = q, ncol = k)
  B <- qr.Q(qr(B))[, 1:k, drop = FALSE]
  
  # Initialize A for convergence check
  A_old <- matrix(0, nrow = p, ncol = k)
  
  for (iter in 1:max_iter) {
    # --- Update A given B ---
    # Solve least squares: A = (X'X)^{-1} X'Y B
    A <- solve(Sxx, crossprod(Xc, Yc) %*% B)
    # Normalize A so that A' Sxx A = I
    norm_A <- sqrt(diag(t(A) %*% Sxx %*% A))
    A <- A %*% diag(1 / norm_A)
    
    # --- Update B given A ---
    # Solve least squares: B = (Y'Y)^{-1} Y'X A
    B <- solve(Syy, crossprod(Yc, Xc) %*% A)
    # Normalize B so that B' Syy B = I
    norm_B <- sqrt(diag(t(B) %*% Syy %*% B))
    B <- B %*% diag(1 / norm_B)
    
    # Check convergence (using the maximum change in A)
    if (max(abs(A - A_old)) < tol) break
    A_old <- A
  }
  
  # Once converged, compute the canonical variates
  U <- Xc %*% A
  V <- Yc %*% B
  
  # Compute the canonical correlations (as the correlation between corresponding variates)
  can_corrs <- sapply(1:k, function(i) cor(U[, i], V[, i]))
  
  return(list(canonical_correlations = can_corrs,
              x_weights = A,
              y_weights = B,
              iterations = iter))
}

# Example usage:
# Suppose X and Y are your data matrices and you want the top 2 canonical directions:
# result <- cca_alternating(X, Y, 2)
# print(result$canonical_correlations)


# Includes lasso penalty
sparse_cca <- function() {
  
}

