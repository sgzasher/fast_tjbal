kmmd_weights <- function(K = NULL, treated_idx, control_idx,
                         eta = 0.001, numdims = NULL,
                         max_numdims = NULL, min_numdims = 1L,
                         newton_tol = 1e-8, newton_maxit = 200L,
                         Phi = NULL) {
  Nco <- length(control_idx)
  Ntr <- length(treated_idx)

  if (!is.null(Phi)) {
    k <- ncol(Phi)
    Phi_co <- Phi[control_idx, , drop = FALSE]
    Phi_tr <- Phi[treated_idx, , drop = FALSE]
    mu_tr <- colMeans(Phi_tr)
    kTc <- as.numeric(Phi_co %*% mu_tr)

    if (Nco == 1L) {
      obj <- sum(Phi_co^2) - 2 * kTc
      return(list(w = 1, objective = obj, converged = TRUE, iterations = 0L, numdims = k))
    }

    sol <- solve_kmmd_dual_phi(Phi_co, kTc, k, eta, newton_tol, newton_maxit)
    sol$numdims <- k
    return(sol)
  }

  if (is.null(K)) stop("Either K or Phi must be provided.")

  if (Nco == 1L) {
    obj <- as.numeric(K[control_idx, control_idx] - 2 * mean(K[treated_idx, control_idx]))
    return(list(w = 1, objective = obj, converged = TRUE, iterations = 0L, numdims = 1L))
  }

  eig <- eigen(K, symmetric = TRUE)
  U   <- eig$vectors
  lam <- eig$values

  pos <- which(lam > 1e-10)
  U   <- U[, pos, drop = FALSE]
  lam <- lam[pos]
  n_pos <- length(pos)

  U_tr <- U[treated_idx, , drop = FALSE]
  U_co <- U[control_idx, , drop = FALSE]
  mu_tr <- colMeans(U_tr)

  if (is.null(max_numdims)) max_numdims <- min(n_pos, Nco - 1L)
  max_numdims <- min(max_numdims, n_pos)

  if (!is.null(numdims)) {
    k_opt <- min(numdims, max_numdims)
    sol <- solve_kmmd_dual_newton(U_co, U_tr, lam, k_opt, eta, newton_tol, newton_maxit)
    sol$numdims <- k_opt
    return(sol)
  }

  best_k   <- min_numdims
  best_bb  <- Inf
  best_sol <- NULL
  rounds_worse <- 0L

  for (k in seq(min_numdims, max_numdims)) {
    sol_k <- solve_kmmd_dual_newton(U_co, U_tr, lam, k, eta, newton_tol, newton_maxit)
    bb <- bias_bound(sol_k$w, U_co, mu_tr, lam, k)

    if (bb < best_bb) {
      best_bb  <- bb
      best_k   <- k
      best_sol <- sol_k
      rounds_worse <- 0L
    } else {
      rounds_worse <- rounds_worse + 1L
    }

    if (bb > 1.25 * best_bb) break
    if (rounds_worse >= 20L) break
  }

  best_sol$numdims    <- best_k
  best_sol$bias_bound <- best_bb
  best_sol
}


#' Solve entropy-regularized KMMD on a rank-k kernel via dual Newton.
#'
#' Primal: min_{w in Delta}  w' K_k w  -  2 kTc' w  +  eta sum(w_i log w_i)
#'   where K_k = U_co_k diag(lam_k) U_co_k'  (rank-k kernel on controls)
#'         kTc_i = sum_j lam_j mu_tr_j U_co_ij  (kernel mean of treated)
#'
#' KKT gives w_i = softmax( (2 kTc - 2 V z) / eta )_i
#'   with V = U_co_k diag(sqrt(lam_k)),  z = V' w
#'
#' Fixed point:  G(z) = z - V' softmax( (2 kTc - 2 V z) / eta ) = 0
#' Newton on this k-dimensional system converges quadratically.
solve_kmmd_dual_newton <- function(U_co, U_tr, lam, k, eta,
                                   tol = 1e-8, maxit = 200L) {
  Nco <- nrow(U_co)
  lam_k  <- lam[1:k]
  U_co_k <- U_co[, 1:k, drop = FALSE]
  mu_tr_k <- colMeans(U_tr[, 1:k, drop = FALSE])

  kTc <- as.numeric(U_co_k %*% (lam_k * mu_tr_k))

  V <- sweep(U_co_k, 2, sqrt(lam_k), `*`)

  z <- rep(0, k)
  converged <- FALSE
  iter <- maxit

  for (i in seq_len(maxit)) {
    Vz <- as.numeric(V %*% z)
    u  <- (2 * kTc - 2 * Vz) / eta
    u  <- u - max(u)
    p  <- exp(u)
    p  <- p / sum(p)

    Gz <- z - as.numeric(crossprod(V, p))
    Gz_norm <- max(abs(Gz))

    if (Gz_norm < tol) {
      converged <- TRUE
      iter <- i
      break
    }

    Vp   <- as.numeric(crossprod(V, p))
    VpVp <- tcrossprod(Vp)
    VdpV <- crossprod(V * sqrt(p))
    J <- diag(k) + (2 / eta) * (VdpV - VpVp)

    dz <- tryCatch(solve(J, Gz), error = function(e) Gz * 0.1)

    alpha <- 1
    for (bt in 1:12) {
      z_try <- z - alpha * dz
      Vz_try <- as.numeric(V %*% z_try)
      u_try <- (2 * kTc - 2 * Vz_try) / eta
      u_try <- u_try - max(u_try)
      p_try <- exp(u_try)
      p_try <- p_try / sum(p_try)
      Gz_try <- z_try - as.numeric(crossprod(V, p_try))

      if (max(abs(Gz_try)) < Gz_norm) break
      alpha <- alpha * 0.5
    }
    z <- z - alpha * dz
  }

  Vz <- as.numeric(V %*% z)
  u  <- (2 * kTc - 2 * Vz) / eta
  u  <- u - max(u)
  w  <- exp(u)
  w  <- w / sum(w)

  Uw <- as.numeric(crossprod(U_co_k, w))
  objective <- as.numeric(sum(lam_k * Uw^2) - 2 * sum(kTc * w) +
    eta * sum(w * log(pmax(w, .Machine$double.xmin))))

  list(w = as.numeric(w), objective = objective, converged = converged, iterations = iter)
}


#' Solve entropy-regularized KMMD directly from feature matrix Phi.
#'
#' V = Phi_co (no eigenvalue scaling needed; K_cc = Phi_co Phi_co').
#' kTc = Phi_co %*% mean(Phi_tr).
solve_kmmd_dual_phi <- function(Phi_co, kTc, k, eta,
                                tol = 1e-8, maxit = 200L) {
  Nco <- nrow(Phi_co)
  V <- Phi_co

  z <- rep(0, k)
  converged <- FALSE
  iter <- maxit

  for (i in seq_len(maxit)) {
    Vz <- as.numeric(V %*% z)
    u  <- (2 * kTc - 2 * Vz) / eta
    u  <- u - max(u)
    p  <- exp(u)
    p  <- p / sum(p)

    Gz <- z - as.numeric(crossprod(V, p))
    Gz_norm <- max(abs(Gz))

    if (Gz_norm < tol) {
      converged <- TRUE
      iter <- i
      break
    }

    Vp   <- as.numeric(crossprod(V, p))
    VpVp <- tcrossprod(Vp)
    VdpV <- crossprod(V * sqrt(p))
    J <- diag(k) + (2 / eta) * (VdpV - VpVp)

    dz <- tryCatch(solve(J, Gz), error = function(e) Gz * 0.1)

    alpha <- 1
    for (bt in 1:12) {
      z_try <- z - alpha * dz
      Vz_try <- as.numeric(V %*% z_try)
      u_try <- (2 * kTc - 2 * Vz_try) / eta
      u_try <- u_try - max(u_try)
      p_try <- exp(u_try)
      p_try <- p_try / sum(p_try)
      Gz_try <- z_try - as.numeric(crossprod(V, p_try))

      if (max(abs(Gz_try)) < Gz_norm) break
      alpha <- alpha * 0.5
    }
    z <- z - alpha * dz
  }

  Vz <- as.numeric(V %*% z)
  u  <- (2 * kTc - 2 * Vz) / eta
  u  <- u - max(u)
  w  <- exp(u)
  w  <- w / sum(w)

  Phiw <- as.numeric(crossprod(Phi_co, w))
  objective <- as.numeric(sum(Phiw^2) - 2 * sum(kTc * w) +
    eta * sum(w * log(pmax(w, .Machine$double.xmin))))

  list(w = as.numeric(w), objective = objective, converged = converged, iterations = iter)
}


bias_bound <- function(w, U_co, mu_tr, lam, k) {
  n_pos <- length(lam)
  if (k >= n_pos) return(0)

  residual_idx <- (k + 1):n_pos
  mu_w_residual <- as.numeric(crossprod(U_co[, residual_idx, drop = FALSE], w))
  imbalance <- mu_tr[residual_idx] - mu_w_residual
  sqrt(sum(lam[residual_idx] * imbalance^2))
}
