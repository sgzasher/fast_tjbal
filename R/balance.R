tjbal_core <- function(
    data,
    Y,
    X,
    Y.match.time = NULL,
    Y.match.npre = NULL,
    Ttot,
    T0,
    id.tr,
    id.co,
    demean = FALSE,
    estimator,
    sigma = NULL,
    eta = 0.001,
    nystrom = NULL,
    nystrom_m = NULL,
    info = TRUE
    ) {

    TT <- length(Ttot)
    Ntr <- length(id.tr)
    Nco <- length(id.co)
    N <- Ntr + Nco
    Tpre <- Ttot[1:T0]
    Tpst <- Ttot[(T0 + 1):TT]

    if (!is.null(Y.match.npre)) {
        if (Y.match.npre == 0) {
            Y.match.time <- NULL
        } else {
            Y.match.time <- Ttot[max(1, (T0 - Y.match.npre + 1)):T0]
        }
    } else {
        if (!is.null(Y.match.time)) {
            Y.match.time <- intersect(Tpre, Y.match.time)
        } else {
            Y.match.time <- Tpre
        }
    }

    data <- data[c(id.tr, id.co), ]
    id.tr <- 1:Ntr
    id.co <- (Ntr + 1):N

    if (demean == TRUE) {
        Y.dm.var <- paste0(Y, ".dm", Ttot)
        Ypre.mean <- apply(data[, paste0(Y, Tpre), drop = FALSE], 1, mean)
        outcome.dm <- data[, paste0(Y, Ttot), drop = FALSE] - matrix(Ypre.mean, N, TT)
        colnames(outcome.dm) <- Y.dm.var
        data <- cbind.data.frame(data, outcome.dm)
        Y.match <- paste0(Y, ".dm", Y.match.time)
        Y.target <- paste0(Y, ".dm", Ttot)
        Y.target.pst <- paste0(Y, ".dm", Tpst)
    } else {
        Y.match <- paste0(Y, Y.match.time)
        Y.target <- paste0(Y, Ttot)
        Y.target.pst <- paste0(Y, Tpst)
    }

    if (is.null(Y.match.time)) {
        Y.match <- NULL
    }
    matchvar <- c(Y.match, X)

    weights.tr <- rep(1 / Ntr, Ntr)
    weights.co <- rep(1 / Nco, Nco)
    w <- c(weights.tr, weights.co * (-1))
    att <- apply(data[, Y.target] * w, 2, sum)
    names(att) <- Ttot
    att.avg <- mean(att[(T0 + 1):TT])

    success <- FALSE

    if (!is.null(matchvar)) {

        if (info == TRUE) {
            cat("Seek balance on:\n")
            cat(paste(matchvar, collapse = ", "), "\n")
            cat("\nOptimization:\n")
        }

        allx <- data.frame(lapply(data[, matchvar, drop = FALSE], as.numeric))
        allx_mat <- as.matrix(allx)

        use_nystrom <- FALSE

        if (estimator == "mean") {
            bal.type <- "mbal"
            K <- linear_kernel_gram(allx_mat)
        } else {
            bal.type <- "kbal"
            use_nystrom <- isTRUE(nystrom) || (is.null(nystrom) && N > 1000)
            if (use_nystrom) {
                if (is.null(nystrom_m)) {
                    nystrom_m <- min(N, ceiling(3 * sqrt(N)))
                }
                if (is.null(sigma)) {
                    dists2 <- as.matrix(dist(allx_mat))^2 / ncol(allx_mat)
                    sigma <- median(sqrt(dists2[dists2 > 0]))
                }
                nys <- kernel_eigen_nystrom(allx_mat, k = nystrom_m, m = nystrom_m, sigma = sigma)
                nys_Phi <- nys$evecs
                K <- NULL
            } else {
                K <- gaussian_kernel_gram(allx_mat, sigma = sigma)
            }
        }

        sol <- tryCatch({
            if (use_nystrom) {
                kmmd_weights(treated_idx = id.tr, control_idx = id.co,
                             eta = eta, Phi = nys_Phi)
            } else {
                kmmd_weights(K, treated_idx = id.tr, control_idx = id.co,
                             eta = eta)
            }
        }, error = function(e) NULL)

        if (!is.null(sol)) {
            success <- TRUE
            weights.co <- sol$w
            weights.tr <- rep(1 / Ntr, Ntr)
            w <- c(weights.tr, weights.co * (-1))
            numdims <- sol$numdims
            if (info == TRUE) {
                conv_str <- if (sol$converged) "converged" else "did not converge"
                cat(paste0("KMMD solver ", conv_str, " in ", sol$iterations,
                    " iterations; objective = ", sprintf("%.6f", sol$objective),
                    if (!is.null(numdims)) paste0("; numdims = ", numdims) else "",
                    "\n"))
            }
        } else {
            success <- FALSE
            message("Solution not found. Equal weights being used.")
        }

    } else {
        bal.type <- "none"
    }

    att <- apply(data[, Y.target] * w, 2, sum)
    names(att) <- Ttot
    att.avg <- mean(att[(T0 + 1):TT])

    Y.var <- paste0(Y, Ttot)
    Y.tr.bar <- apply(data[id.tr, Y.var, drop = FALSE], 2, mean, na.rm = TRUE)
    Y.co.bar <- apply(data[id.co, Y.var, drop = FALSE], 2, mean, na.rm = TRUE)
    Y.ct.bar <- Y.tr.bar - att
    Y.bar <- cbind(Y.tr.bar, Y.ct.bar, Y.co.bar)

    out <- list(
        data.wide = data,
        w = w,
        weights.co = weights.co,
        matchvar = matchvar,
        Y.var = Y.var,
        Y.target = Y.target,
        Y.target.pst = Y.target.pst,
        Y.bar = Y.bar,
        att = att,
        att.avg = att.avg,
        ntreated = rep(Ntr, TT),
        bal.type = bal.type
    )

    if (!is.null(matchvar)) {
        out$success <- success
        if (success == TRUE) {
            bal.table <- compute_balance_table(data, matchvar, id.tr, id.co, weights.co)
            out$K <- K
            out$sigma <- if (!is.null(K)) attr(K, "sigma") else sigma
            out$numdims <- numdims
            out$bal.table <- bal.table
        }
    }

    return(out)
}


compute_balance_table <- function(data, matchvar, id.tr, id.co, weights.co) {

    Ntr <- length(id.tr)

    weighted.sd <- function(vec, w) {
        sqrt(sum(w * (vec - weighted.mean(vec, w))^2))
    }

    if (Ntr > 1) {
        mean.tr <- apply(data[id.tr, matchvar, drop = FALSE], 2, mean)
        sd.tr <- apply(data[id.tr, matchvar, drop = FALSE], 2, sd)
        mean.co.pre <- apply(data[id.co, matchvar, drop = FALSE], 2, mean)
        sd.co.pre <- apply(data[id.co, matchvar, drop = FALSE], 2, sd)
        mean.co.pst <- apply(data[id.co, matchvar, drop = FALSE], 2, weighted.mean, weights.co)
        sd.co.pst <- apply(data[id.co, matchvar, drop = FALSE], 2, weighted.sd, weights.co)
        diff.pre <- (mean.tr - mean.co.pre) / sd.tr
        diff.pst <- (mean.tr - mean.co.pst) / sd.tr
        bal.table <- cbind.data.frame(
            mean.tr, mean.co.pre, mean.co.pst,
            sd.tr, sd.co.pre, sd.co.pst,
            diff.pre, diff.pst
        )
    } else {
        mean.tr <- apply(data[id.tr, matchvar, drop = FALSE], 2, mean)
        mean.co.pre <- apply(data[id.co, matchvar, drop = FALSE], 2, mean)
        mean.co.pst <- apply(data[id.co, matchvar, drop = FALSE], 2, weighted.mean, weights.co)
        diff.pre <- (mean.tr - mean.co.pre) / abs(mean.tr)
        diff.pst <- (mean.tr - mean.co.pst) / abs(mean.tr)
        bal.table <- cbind.data.frame(
            mean.tr, mean.co.pre, mean.co.pst,
            diff.pre, diff.pst
        )
    }

    return(bal.table)
}
