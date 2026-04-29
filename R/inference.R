jackknife_att <- function(data, Y.target, id.tr, id.co, K, T0, TT, eta, nsims, parallel, cores) {
  Ntr <- length(id.tr)
  Nco <- length(id.co)
  N <- Ntr + Nco
  njacks <- min(nsims, Ntr)
  drop.id <- sample(id.tr, njacks, replace = FALSE)

  one.jack <- function(drop) {
    sample.id <- c(id.co, setdiff(id.tr, drop))
    Ntr_jack <- Ntr - 1
    N_jack <- Ntr_jack + Nco
    w.jack <- rep(1 / Ntr_jack, N_jack)

    if (!is.null(K)) {
      sol <- kmmd_weights(K[sample.id, sample.id, drop = FALSE],
                          treated_idx = (Nco + 1):N_jack,
                          control_idx = 1:Nco,
                          eta = eta)
      w.jack[1:Nco] <- sol$w * (-1)
    } else {
      w.jack[1:Nco] <- rep(-1 / Nco, Nco)
    }

    att <- apply(data[sample.id, Y.target] * w.jack, 2, sum)
    att.avg <- mean(att[(T0 + 1):TT])
    list(att = att, att.avg = att.avg)
  }

  att.sims <- matrix(0, TT, njacks)
  att.avg.sims <- matrix(0, njacks, 1)

  if (parallel && njacks >= 8) {
    if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    results <- foreach::foreach(j = 1:njacks, .inorder = FALSE,
                                .export = c("kmmd_weights", "solve_kmmd_dual_newton", "solve_kmmd_dual_phi", "bias_bound")) %dopar% {
      one.jack(drop.id[j])
    }
    parallel::stopCluster(cl)
    for (j in seq_len(njacks)) {
      att.sims[, j] <- results[[j]]$att
      att.avg.sims[j, ] <- results[[j]]$att.avg
    }
  } else {
    for (j in seq_len(njacks)) {
      res <- one.jack(drop.id[j])
      att.sims[, j] <- res$att
      att.avg.sims[j, ] <- res$att.avg
    }
  }

  list(att.sims = att.sims, att.avg.sims = att.avg.sims, njacks = njacks)
}


bootstrap_att <- function(data, Y.target, id.tr, id.co, K, T0, TT, eta, nsims, parallel, cores) {
  Ntr <- length(id.tr)
  Nco <- length(id.co)
  N <- Ntr + Nco

  one.boot <- function() {
    sample.id <- c(sample(id.co, Nco, replace = TRUE), sample(id.tr, Ntr, replace = TRUE))
    w.boot <- rep(1 / Ntr, N)

    if (!is.null(K)) {
      sol <- kmmd_weights(K[sample.id, sample.id, drop = FALSE],
                          treated_idx = (Nco + 1):N,
                          control_idx = 1:Nco,
                          eta = eta)
      w.boot[1:Nco] <- sol$w * (-1)
    } else {
      w.boot[1:Nco] <- rep(-1 / Nco, Nco)
    }

    att <- apply(data[sample.id, Y.target] * w.boot, 2, sum)
    att.avg <- mean(att[(T0 + 1):TT])
    list(att = att, att.avg = att.avg)
  }

  att.sims <- matrix(0, TT, nsims)
  att.avg.sims <- matrix(0, nsims, 1)

  if (parallel && nsims >= 8) {
    if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    results <- foreach::foreach(j = 1:nsims, .inorder = FALSE,
                                .export = c("kmmd_weights", "solve_kmmd_dual_newton", "solve_kmmd_dual_phi", "bias_bound")) %dopar% {
      one.boot()
    }
    parallel::stopCluster(cl)
    for (j in seq_len(nsims)) {
      att.sims[, j] <- results[[j]]$att
      att.avg.sims[j, ] <- results[[j]]$att.avg
    }
  } else {
    for (j in seq_len(nsims)) {
      res <- one.boot()
      att.sims[, j] <- res$att
      att.avg.sims[j, ] <- res$att.avg
    }
  }

  list(att.sims = att.sims, att.avg.sims = att.avg.sims)
}


fixed_weight_att <- function(data, Y.target, w, id.tr, id.co, T0, TT, nsims) {
  Ntr <- length(id.tr)
  Nco <- length(id.co)
  N <- Ntr + Nco
  att.sims <- matrix(0, TT, nsims)
  att.avg.sims <- matrix(0, nsims, 1)

  for (j in 1:nsims) {
    sample.id <- c(sample(id.tr, Ntr, replace = TRUE),
                   sample(id.co, Nco, replace = TRUE))
    w.boot <- w[sample.id]
    w.boot[1:Ntr] <- w.boot[1:Ntr] / sum(w.boot[1:Ntr])
    w.boot[(Ntr + 1):N] <- w.boot[(Ntr + 1):N] / sum(w.boot[(Ntr + 1):N]) * (-1)
    att.sims[, j] <- apply(data[sample.id, Y.target] * w.boot, 2, sum)
    att.avg.sims[j, ] <- mean(att.sims[(T0 + 1):TT, j])
  }

  list(att.sims = att.sims, att.avg.sims = att.avg.sims)
}


compute_se_ci <- function(att, att.sims, att.avg, att.avg.sims, Ntr, TT, conf.lvl, vce_type) {
  se.att <- apply(att.sims, 1, function(vec) sd(vec, na.rm = TRUE))
  se.att.avg <- sd(att.avg.sims, na.rm = TRUE)

  if (vce_type == "jackknife") {
    se.att <- se.att * (Ntr - 1) / sqrt(Ntr)
    se.att.avg <- se.att.avg * (Ntr - 1) / sqrt(Ntr)
  }

  z.att <- att / se.att
  z.att.avg <- att.avg / se.att.avg

  pvalue.att <- (1 - pnorm(abs(z.att))) * 2
  pvalue.att.avg <- (1 - pnorm(abs(z.att.avg))) * 2

  c.value <- qnorm(0.5 + conf.lvl / 2)
  CI.att <- cbind(att - c.value * se.att, att + c.value * se.att)
  CI.att.avg <- c(att.avg - c.value * se.att.avg, att.avg + c.value * se.att.avg)

  est.att <- cbind(att, se.att, z.att, CI.att, pvalue.att, ntreated = rep(Ntr, TT))
  est.att[abs(est.att) < 1e-5] <- 0
  colnames(est.att) <- c("ATT", "S.E.", "z-score", "CI.lower", "CI.upper", "p.value", "n.Treated")

  est.att.avg <- t(as.matrix(c(att.avg, se.att.avg, z.att.avg, CI.att.avg, pvalue.att.avg)))
  colnames(est.att.avg) <- c("ATT", "S.E.", "z-score", "CI.lower", "CI.upper", "p.value")

  list(
    est.att = round(est.att, 4),
    est.att.avg = round(est.att.avg, 4),
    att.sims = att.sims,
    att.avg.sims = att.avg.sims
  )
}
