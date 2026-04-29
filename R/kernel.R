gaussian_kernel_gram <- function(X, sigma = NULL, normalize_n = FALSE) {
  X <- as.matrix(X)
  n <- nrow(X)

  dists2 <- as.matrix(dist(X))^2 / ncol(X)

  if (is.null(sigma)) {
    dists <- sqrt(dists2)
    sigma <- median(dists[dists > 0])
  }

  G <- exp(-dists2 / (2 * sigma^2))
  attr(G, "sigma") <- sigma

  if (normalize_n) G / n else G
}

gaussian_kernel_cross <- function(X, Y, sigma) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  X_norm2 <- rowSums(X^2)
  Y_norm2 <- rowSums(Y^2)
  cross   <- X %*% t(Y)
  dists2  <- (outer(X_norm2, Y_norm2, "+") - 2 * cross) / ncol(X)

  exp(-dists2 / (2 * sigma^2))
}

linear_kernel_gram <- function(X) {
  X <- as.matrix(X)
  X %*% t(X)
}

kernel_eigen_exact <- function(X, k, sigma = NULL, ridge_mult = 0, tol = 1e-10) {
  X <- as.matrix(X)
  n <- nrow(X)

  K <- gaussian_kernel_gram(X, sigma = sigma, normalize_n = TRUE)
  sigma <- attr(K, "sigma")

  if (ridge_mult > 0) {
    diag(K) <- diag(K) + ridge_mult / n
  }

  K <- 0.5 * (K + t(K))
  eig  <- eigen(K, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors

  pos_idx <- which(vals > tol)
  if (length(pos_idx) == 0L) stop("No positive eigenvalues found in K.")
  k_eff <- min(as.integer(k), length(pos_idx))
  idx   <- pos_idx[seq_len(k_eff)]

  list(
    evecs = vecs[, idx, drop = FALSE],
    evals = vals[idx],
    sigma = sigma
  )
}

kernel_eigen_nystrom <- function(X, k, m, sigma = NULL, landmark_idx = NULL, jitter = 1e-8) {
  X <- as.matrix(X)
  n <- nrow(X)

  if (is.null(landmark_idx)) {
    m <- min(m, n)
    landmark_idx <- sample.int(n, size = m, replace = FALSE)
  } else {
    landmark_idx <- unique(landmark_idx)
    m <- length(landmark_idx)
  }

  X_land <- X[landmark_idx, , drop = FALSE]

  if (is.null(sigma)) {
    dists_land <- sqrt(as.matrix(dist(X_land))^2 / ncol(X_land))
    sigma <- median(dists_land[dists_land > 0])
  }

  dists2_mm <- as.matrix(dist(X_land))^2 / ncol(X_land)
  K_mm <- exp(-dists2_mm / (2 * sigma^2))
  attr(K_mm, "sigma") <- sigma
  K_mm <- 0.5 * (K_mm + t(K_mm))
  diag(K_mm) <- diag(K_mm) + jitter

  X_norm2     <- rowSums(X^2)
  Land_norm2  <- rowSums(X_land^2)
  cross       <- X %*% t(X_land)
  dists2_nm   <- (outer(X_norm2, Land_norm2, "+") - 2 * cross) / ncol(X)
  C <- exp(-dists2_nm / (2 * sigma^2))

  eig   <- eigen(K_mm, symmetric = TRUE)
  lam_m <- eig$values
  U_m   <- eig$vectors

  pos_idx <- which(lam_m > 1e-10)
  if (length(pos_idx) == 0L) {
    stop("No positive eigenvalues found in K_mm (after jitter); sigma might be pathological.")
  }

  k_eff <- min(k, length(pos_idx), m)
  idx   <- pos_idx[seq_len(k_eff)]

  lam_k        <- lam_m[idx]
  U_k          <- U_m[, idx, drop = FALSE]
  lam_k_reg    <- pmax(lam_k, 1e-8)
  Lambda_inv_sqrt <- diag(1 / sqrt(lam_k_reg), nrow = k_eff, ncol = k_eff)

  Phi <- C %*% U_k %*% Lambda_inv_sqrt
  colnames(Phi) <- paste0("phi", seq_len(ncol(Phi)))

  list(
    evecs        = Phi,
    evals        = lam_k / n,
    sigma        = sigma,
    landmark_idx = landmark_idx
  )
}
