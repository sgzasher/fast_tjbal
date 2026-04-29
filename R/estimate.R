tjbal_estimate_single <- function(
    data,
    Y,
    D,
    X,
    Y.match.time = NULL,
    Y.match.npre = NULL,
    Ttot,
    demean = FALSE,
    estimator = "kernel",
    sigma = NULL,
    eta = 0.001,
    nystrom = NULL,
    nystrom_m = NULL,
    print.baltable = FALSE,
    vce = "jackknife",
    conf.lvl = 0.95,
    nsims = 500,
    parallel = TRUE,
    cores = NULL
    ) {

    TT <- length(Ttot)
    id.tr <- which(data$treat == 1)
    id.co <- which(data$treat == 0)
    Ntr <- length(id.tr)
    Nco <- length(id.co)
    N <- Ntr + Nco
    T0 <- unique(data$T0[id.tr])
    Tpre <- Ttot[1:T0]
    Tpst <- Ttot[(T0 + 1):TT]

    out <- tjbal_core(
        data = data, Y = Y, X = X, Y.match.time = Y.match.time,
        Y.match.npre = Y.match.npre, Ttot = Ttot, T0 = T0,
        id.tr = id.tr, id.co = id.co,
        demean = demean, estimator = estimator, sigma = sigma,
        eta = eta, nystrom = nystrom, nystrom_m = nystrom_m,
        info = TRUE
    )

    att <- out$att
    att.avg <- out$att.avg
    Y.var <- out$Y.var
    Y.bar <- out$Y.bar
    w <- out$w
    weights.co <- out$weights.co
    matchvar <- out$matchvar
    Y.target <- out$Y.target
    Y.target.pst <- out$Y.target.pst
    K <- out$K
    data <- out$data.wide
    id.tr <- 1:Ntr
    id.co <- (Ntr + 1):N

    if (print.baltable == TRUE && !is.null(matchvar) && !is.null(out$bal.table)) {
        cat("\nBalance Table\n")
        print(round(out$bal.table, 4))
    }

    if (vce == "jackknife") {
        cat("\nJackknife...\n")
        inf <- jackknife_att(
            data = data, Y.target = Y.target,
            id.tr = id.tr, id.co = id.co,
            K = K, T0 = T0, TT = TT, eta = eta,
            nsims = nsims, parallel = parallel, cores = cores
        )
        se_out <- compute_se_ci(
            att, inf$att.sims, att.avg, inf$att.avg.sims,
            Ntr, TT, conf.lvl, "jackknife"
        )
        out <- c(out, se_out)
    } else if (vce == "bootstrap") {
        cat("\nBootstrapping...\n")
        inf <- bootstrap_att(
            data = data, Y.target = Y.target,
            id.tr = id.tr, id.co = id.co,
            K = K, T0 = T0, TT = TT, eta = eta,
            nsims = nsims, parallel = parallel, cores = cores
        )
        se_out <- compute_se_ci(
            att, inf$att.sims, att.avg, inf$att.avg.sims,
            Ntr, TT, conf.lvl, "bootstrap"
        )
        out <- c(out, se_out)
    } else if (vce == "fixed.weights") {
        inf <- fixed_weight_att(
            data = data, Y.target = Y.target, w = w,
            id.tr = id.tr, id.co = id.co,
            T0 = T0, TT = TT, nsims = nsims
        )
        se_out <- compute_se_ci(
            att, inf$att.sims, att.avg, inf$att.avg.sims,
            Ntr, TT, conf.lvl, "fixed.weights"
        )
        out <- c(out, se_out)
    }

    out <- c(out, list(
        id.tr = id.tr,
        id.co = id.co,
        Y.var = Y.var,
        Ttot = Ttot,
        Tpre = Tpre,
        Tpst = Tpst,
        T0 = T0,
        N = N,
        Ntr = Ntr,
        Nco = Nco,
        ntreated = rep(Ntr, TT)
    ))

    return(out)
}


tjbal_estimate_multi <- function(
    data,
    Y,
    D,
    X,
    Y.match.time = NULL,
    Y.match.npre = NULL,
    Ttot,
    demean = FALSE,
    estimator = "kernel",
    sigma = NULL,
    eta = 0.001,
    nystrom = NULL,
    nystrom_m = NULL,
    vce = "jackknife",
    conf.lvl = 0.95,
    nsims = 200,
    parallel = TRUE,
    cores = NULL
    ) {

    TT <- length(Ttot)
    id.tr <- which(data$treat == 1)
    id.co <- which(data$treat == 0)
    Ntr <- length(id.tr)
    Nco <- length(id.co)
    N <- Ntr + Nco
    Y.var <- paste0(Y, Ttot)
    units <- as.character(unique(data$unit))

    T0.all <- data$T0
    names(T0.all) <- units
    T0.tr <- T0.all[id.tr]
    tb.T0 <- table(T0.tr)
    T0.unique <- as.numeric(names(tb.T0))
    T0.count <- as.numeric(tb.T0)
    T0.max <- max(T0.tr)
    T0.min <- min(T0.tr)
    T0.names <- paste0("T0 = ", Ttot[T0.unique])
    nT0 <- length(T0.unique)

    time.adj <- c(-(T0.max - 1):(TT - T0.min))
    time.adj.max <- TT - T0.min
    TT.adj <- length(time.adj)

    sub.weights.co <- matrix(NA, length(id.co), nT0)
    rownames(sub.weights.co) <- units[id.co]
    colnames(sub.weights.co) <- T0.names
    sub.att <- sub.Ytr.avg <- sub.Yct.avg <- matrix(NA, TT, nT0)
    sub.att.avg <- rep(NA, nT0)
    rownames(sub.att) <- rownames(sub.Ytr.avg) <- rownames(sub.Yct.avg) <- Ttot
    names(sub.att.avg) <- colnames(sub.att) <- colnames(sub.Ytr.avg) <- colnames(sub.Yct.avg) <- T0.names
    sub.Ytr.adj <- sub.ntr.pst <- sub.ntr <- sub.att.adj <- matrix(0, TT.adj, nT0)
    se.sub.att.adj <- matrix(NA, TT.adj, nT0)
    rownames(sub.Ytr.adj) <- rownames(sub.ntr.pst) <- rownames(sub.ntr) <- rownames(se.sub.att.adj) <- rownames(sub.att.adj) <- time.adj
    colnames(sub.Ytr.adj) <- colnames(sub.ntr.pst) <- colnames(sub.ntr) <- colnames(se.sub.att.adj) <- colnames(sub.att.adj) <- T0.names

    K.list <- data.list <- matchvar.list <- bal.table.list <- Y.target.list <- vector("list", length = nT0)
    success <- rep(0, nT0)
    bias.ratios <- rep(1, nT0)

    cat("Balancing...\n")
    for (i in 1:nT0) {
        T0_i <- T0.unique[i]
        cat(paste0("Subgroup T0 = ", T0_i, ": "))
        id.tr.one <- which(T0.all == T0_i)

        core.out <- tjbal_core(
            data = data, Y = Y, X = X, Y.match.time = Y.match.time,
            Y.match.npre = Y.match.npre, Ttot = Ttot, T0 = T0_i,
            id.tr = id.tr.one, id.co = id.co,
            demean = demean, estimator = estimator, sigma = sigma,
            eta = eta, nystrom = nystrom, nystrom_m = nystrom_m,
            info = FALSE
        )

        att_i <- core.out$att
        data.list[[i]] <- core.out$data.wide
        Y.target.list[[i]] <- core.out$Y.target

        if (!is.null(core.out$matchvar)) {
            matchvar.list[[i]] <- core.out$matchvar
            success[i] <- core.out$success
            if (core.out$success) {
                K.list[[i]] <- core.out$K
            }
        }

        sub.weights.co[, i] <- core.out$weights.co
        sub.att[, i] <- att_i
        sub.att.avg[i] <- core.out$att.avg
        sub.Ytr.avg[, i] <- core.out$Y.bar[, 1]
        sub.Yct.avg[, i] <- core.out$Y.bar[, 2]

        fill.start <- T0.max - T0_i + 1
        fill.end <- fill.start + length(att_i) - 1
        sub.ntr[fill.start:fill.end, i] <- T0.count[i]
        sub.att.adj[fill.start:fill.end, i] <- att_i
        sub.Ytr.adj[fill.start:fill.end, i] <- sub.Ytr.avg[, i]

        if (!is.null(core.out$matchvar)) {
            bal.table.list[[i]] <- core.out$bal.table
        }
    }

    ntreated <- rowSums(sub.ntr)
    att <- rowSums(sub.att.adj * sub.ntr) / ntreated
    names(ntreated) <- names(att) <- time.adj

    Y.bar.tr <- rowSums(sub.Ytr.adj * sub.ntr) / ntreated
    Y.bar.ct <- Y.bar.tr - att
    Y.bar <- cbind(Y.bar.tr, Y.bar.ct)

    sub.ntr.pst <- sub.ntr
    sub.ntr.pst[time.adj <= 0, ] <- 0

    att.avg <- sum(sub.att.adj * sub.ntr.pst, na.rm = TRUE) / sum(sub.ntr.pst)

    weights.co <- apply(sub.weights.co, 1, weighted.mean, T0.count)
    names(weights.co) <- units[id.co]

    group.stats <- cbind(
        T0 = T0.unique, time = Ttot[T0.unique],
        Ntr = T0.count, success = success
    )

    if (vce == "bootstrap") {
        stop("Bootstrap is not supported for staggered adoption designs. Use vce = 'jackknife'.")
    }

    if (vce == "jackknife") {
        njacks <- min(nsims, Ntr)
        cat("\nJackknife...\n")

        drop.id.pos <- sample(1:Ntr, njacks, replace = FALSE)
        drop.id <- id.tr[drop.id.pos]
        drop.id.T0s <- T0.tr[drop.id.pos]
        drop.id.list <- split(drop.id, drop.id.T0s)

        sub.att.jack <- vector("list", length = nT0)
        sub.att.avg.jack <- vector("list", length = nT0)
        names(sub.att.avg.jack) <- names(sub.att.jack) <- colnames(sub.att)

        if (parallel) {
            if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)
            para.clusters <- parallel::makeCluster(cores)
            doParallel::registerDoParallel(para.clusters)
        }

        for (i in 1:nT0) {
            T0_i <- T0.unique[i]
            cat("\nDropping units from Subgroup T0 =", T0_i)
            drop.id.oneT0 <- drop.id.list[[as.character(T0_i)]]

            if (is.null(drop.id.oneT0) || length(drop.id.oneT0) == 0) next

            nid <- length(drop.id.oneT0)
            Ntr_sub <- T0.count[i]
            data_sub <- data.list[[i]]
            K_sub <- K.list[[i]]
            Y.target_sub <- Y.target.list[[i]]
            N_sub <- nrow(data_sub)

            sub.att.oneT0 <- matrix(NA, TT, nid)
            sub.att.avg.oneT0 <- rep(NA, nid)

            if (Ntr_sub <= 1) {
                sub.att.jack[[i]] <- sub.att.oneT0
                sub.att.avg.jack[[i]] <- sub.att.avg.oneT0
                next
            }

            one.jack.sub <- function(orig_id) {
                row_drop <- which(data_sub$id == orig_id)
                sample.id <- setdiff(1:N_sub, row_drop)
                N_jack <- length(sample.id)
                Ntr_jack <- Ntr_sub - 1

                w.jack <- rep(1 / Ntr_jack, N_jack)

                if (!is.null(K_sub)) {
                    tr_in_sample <- setdiff(1:Ntr_sub, row_drop)
                    co_in_sample <- (Ntr_sub + 1):N_sub
                    new_tr_pos <- match(tr_in_sample, sample.id)
                    new_co_pos <- match(co_in_sample, sample.id)

                    sol <- tryCatch(
                        kmmd_weights(
                            K_sub[sample.id, sample.id, drop = FALSE],
                            treated_idx = new_tr_pos,
                            control_idx = new_co_pos,
                            eta = eta
                        ),
                        error = function(e) NULL
                    )
                    if (!is.null(sol)) {
                        w.jack[new_co_pos] <- sol$w * (-1)
                    } else {
                        w.jack[new_co_pos] <- rep(-1 / Nco, Nco)
                    }
                } else {
                    w.jack[(Ntr_jack + 1):N_jack] <- rep(-1 / Nco, Nco)
                }

                att <- apply(data_sub[sample.id, Y.target_sub] * w.jack, 2, sum)
                att.avg <- mean(att[(T0_i + 1):TT])
                list(att = att, att.avg = att.avg)
            }

            if (parallel && nid >= 8) {
                jack.out <- foreach::foreach(
                    j = 1:nid, .inorder = FALSE,
                    .export = c("kmmd_weights", "solve_kmmd_dual_newton", "solve_kmmd_dual_phi", "bias_bound")
                ) %dopar% {
                    one.jack.sub(drop.id.oneT0[j])
                }
                for (j in 1:nid) {
                    sub.att.oneT0[, j] <- jack.out[[j]]$att
                    sub.att.avg.oneT0[j] <- jack.out[[j]]$att.avg
                }
            } else {
                for (j in 1:nid) {
                    res <- one.jack.sub(drop.id.oneT0[j])
                    sub.att.oneT0[, j] <- res$att
                    sub.att.avg.oneT0[j] <- res$att.avg
                    if (j %% 50 == 0) cat(".")
                }
            }

            sub.att.jack[[i]] <- sub.att.oneT0
            sub.att.avg.jack[[i]] <- sub.att.avg.oneT0
        }

        if (parallel) parallel::stopCluster(para.clusters)

        c.value <- qnorm(0.5 + conf.lvl / 2)

        CI.sub.att.low <- CI.sub.att.high <- pvalue.sub.att <- z.sub.att <- se.sub.att <- matrix(NA, TT, nT0)
        se.sub.att.avg <- rep(NA, nT0)
        rownames(se.sub.att) <- Ttot
        names(se.sub.att.avg) <- colnames(se.sub.att) <- T0.names

        for (i in 1:nT0) {
            if (is.null(sub.att.jack[[i]])) next
            Ntr.oneT0 <- ncol(sub.att.jack[[i]])
            if (Ntr.oneT0 > 1) {
                se.sub.att[, i] <- apply(sub.att.jack[[i]], 1, function(vec) sd(vec, na.rm = TRUE)) * (Ntr.oneT0 - 1) / sqrt(Ntr.oneT0)
                se.sub.att.avg[i] <- sd(sub.att.avg.jack[[i]], na.rm = TRUE) * (Ntr.oneT0 - 1) / sqrt(Ntr.oneT0)
                z.sub.att[, i] <- sub.att[, i] / se.sub.att[, i]
                pvalue.sub.att[, i] <- (1 - pnorm(abs(z.sub.att[, i]))) * 2
                CI.sub.att.low[, i] <- sub.att[, i] - c.value * se.sub.att[, i]
                CI.sub.att.high[, i] <- sub.att[, i] + c.value * se.sub.att[, i]
            }
        }

        group.above1 <- which(!is.na(se.sub.att.avg))
        for (i in group.above1) {
            T0_i <- T0.unique[i]
            fill.start <- T0.max - T0_i + 1
            fill.end <- fill.start + TT - 1
            se.sub.att.adj[fill.start:fill.end, i] <- se.sub.att[, i]
        }

        vce.sub.att.adj <- se.sub.att.adj^2
        sum.vce.bytime <- apply(
            vce.sub.att.adj[, group.above1, drop = FALSE] * sub.ntr[, group.above1, drop = FALSE]^2,
            1, sum, na.rm = TRUE
        )
        sum.obs.bytime <- apply(sub.ntr[, group.above1, drop = FALSE], 1, sum)
        sum.vce.bytime[sum.vce.bytime == 0] <- NA
        sum.obs.bytime[sum.obs.bytime == 0] <- NA
        se.att <- sqrt(sum.vce.bytime / sum.obs.bytime^2)

        sub.obs <- apply(sub.ntr.pst[, group.above1, drop = FALSE], 2, sum)
        vce.sub.att.avg <- se.sub.att.avg[group.above1]^2
        se.att.avg <- sqrt(sum(vce.sub.att.avg * sub.obs^2) / (sum(sub.obs)^2))

        z.att <- att / se.att
        z.att.avg <- att.avg / se.att.avg
        z.sub.att.avg <- sub.att.avg / se.sub.att.avg

        pvalue.att <- (1 - pnorm(abs(z.att))) * 2
        pvalue.att.avg <- (1 - pnorm(abs(z.att.avg))) * 2
        pvalue.sub.att.avg <- (1 - pnorm(abs(z.sub.att.avg))) * 2

        CI.att <- cbind(att - c.value * se.att, att + c.value * se.att)
        CI.att.avg <- c(att.avg - c.value * se.att.avg, att.avg + c.value * se.att.avg)
        CI.sub.att.avg <- cbind(
            sub.att.avg - c.value * se.sub.att.avg,
            sub.att.avg + c.value * se.sub.att.avg
        )

        est.names <- c("ATT", "S.E.", "z-score", "CI.lower", "CI.upper", "p.value", "n.Treated")

        est.att <- cbind(att, se.att, z.att, CI.att, pvalue.att, ntreated = ntreated)
        est.att[abs(est.att) < 1e-5] <- 0
        colnames(est.att) <- est.names

        est.att.avg <- t(as.matrix(c(att.avg, se.att.avg, z.att.avg, CI.att.avg, pvalue.att.avg)))
        colnames(est.att.avg) <- est.names[1:6]

        est.sub.att.avg <- cbind(
            sub.att.avg, se.sub.att.avg, z.sub.att.avg,
            CI.sub.att.avg, pvalue.sub.att.avg, ntreated = T0.count
        )
        est.sub.att.avg[abs(est.sub.att.avg) < 1e-5] <- 0
        colnames(est.sub.att.avg) <- est.names

        est.sub.att <- array(NA, dim = c(TT, 7, nT0),
            dimnames = list(Ttot, est.names, paste("T0 =", Ttot[T0.unique])))
        for (i in 1:nT0) {
            est.sub.att[,, i] <- cbind(
                sub.att[, i], se.sub.att[, i], z.sub.att[, i],
                CI.sub.att.low[, i], CI.sub.att.high[, i],
                pvalue.sub.att[, i], rep(T0.count[i], TT)
            )
        }

        out.inference <- list(
            est.att = round(est.att, 4),
            est.att.avg = round(est.att.avg, 4),
            est.sub.att = round(est.sub.att, 4),
            est.sub.att.avg = round(est.sub.att.avg, 4)
        )
    }

    sub.att.adj[sub.ntr == 0] <- NA

    out <- list(
        data.wide = data,
        id.tr = id.tr,
        id.co = id.co,
        Y.var = Y.var,
        matchvar.list = matchvar.list,
        Ttot = Ttot,
        N = N,
        Ntr = Ntr,
        Nco = Nco,
        T0 = T0.unique,
        T0.all = T0.all,
        T0.tr = T0.tr,
        weights.co = weights.co,
        Y.bar = Y.bar,
        sub.weights.co = sub.weights.co,
        sub.Ytr.avg = sub.Ytr.avg,
        sub.Yct.avg = sub.Yct.avg,
        sub.att = sub.att,
        sub.ntr = sub.ntr,
        sub.att.adj = sub.att.adj,
        ntreated = ntreated,
        att = att,
        att.avg = att.avg,
        success = success,
        group.stats = group.stats,
        bal.table.list = bal.table.list
    )

    if (vce == "jackknife") {
        out <- c(out, out.inference)
    }

    return(out)
}
