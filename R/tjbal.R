tjbal <- function(...) {
    UseMethod("tjbal")
}

tjbal.formula <- function(
    formula = NULL,
    data,
    X.avg.time = NULL,
    index,
    trim.npre = 0,
    Y.match.time = NULL,
    Y.match.npre = NULL,
    demean = TRUE,
    estimator = "kernel",
    sigma = NULL,
    eta = 0.001,
    nystrom = NULL,
    nystrom_m = NULL,
    print.baltable = TRUE,
    vce = "jackknife",
    conf.lvl = 0.95,
    nsims = 200,
    parallel = TRUE,
    cores = NULL,
    seed = NULL,
    ...
    ) {
    varnames <- all.vars(formula)
    Yname <- varnames[1]
    Dname <- varnames[2]
    Xname <- if (length(varnames) > 2) varnames[3:length(varnames)] else NULL

    namesData <- colnames(data)
    for (i in 1:length(varnames)) {
        if (!varnames[i] %in% namesData) {
            stop(paste0("variable \"", varnames[i], "\" is not in the dataset."))
        }
    }

    out <- tjbal.default(data = data, Y = Yname, D = Dname, X = Xname,
        X.avg.time = X.avg.time, index = index, trim.npre = trim.npre,
        Y.match.time = Y.match.time, Y.match.npre = Y.match.npre,
        demean = demean, estimator = estimator, sigma = sigma,
        eta = eta, nystrom = nystrom, nystrom_m = nystrom_m,
        print.baltable = print.baltable,
        vce = vce, conf.lvl = conf.lvl, nsims = nsims,
        parallel = parallel, cores = cores, seed = seed)

    out$call <- match.call()
    out$formula <- formula
    return(out)
}

tjbal.default <- function(
    data,
    Y,
    D,
    X = NULL,
    X.avg.time = NULL,
    index,
    trim.npre = 0,
    Y.match.time = NULL,
    Y.match.npre = NULL,
    demean = TRUE,
    estimator = "kernel",
    sigma = NULL,
    eta = 0.001,
    nystrom = NULL,
    nystrom_m = NULL,
    print.baltable = TRUE,
    vce = "jackknife",
    conf.lvl = 0.95,
    nsims = 200,
    parallel = TRUE,
    cores = NULL,
    seed = NULL,
    ...
    ) {

    if (!is.null(seed)) set.seed(seed)
    if (inherits(data, "tbl_df")) data <- as.data.frame(data)
    if (!is.data.frame(data)) stop("Not a data frame")
    data <- droplevels(data)

    if (length(index) != 2 || sum(index %in% colnames(data)) != 2) {
        stop("\"index\" option misspecified. Try, for example, index = c(\"unit.id\", \"time\").")
    }

    if (vce == "boot") vce <- "bootstrap"
    if (vce == "jack") vce <- "jackknife"
    if (vce == "fixed") vce <- "fixed.weights"

    if (!estimator %in% c("mean", "kernel")) {
        stop("\"estimator\" must be either \"mean\" or \"kernel\".")
    }

    if (!is.null(Y.match.time)) {
        if (Y.match.time[1] == "none") {
            Y.match.time <- NULL
        }
    }

    if (!is.null(sigma) && !is.numeric(sigma)) {
        stop("\"sigma\" needs to be numeric; the default is 2.")
    }

    if (is.null(nsims)) nsims <- 200

    Yname <- Y
    Dname <- D
    Xname <- X

    id <- index[1]
    time <- index[2]
    TT <- length(unique(data[, time]))
    N <- length(unique(data[, id]))
    p <- length(Xname)

    if (var(table(data[, id])) + var(table(data[, time])) > 0) {
        stop("The panel is not balanced.")
    }
    if (!is.numeric(data[, time])) {
        stop("The time indicator must be numeric.")
    }

    if (sum(is.na(data[, Yname])) > 0) stop(paste("Missing values in variable \"", Yname, "\".", sep = ""))
    if (sum(is.na(data[, Dname])) > 0) stop(paste("Missing values in variable \"", Dname, "\".", sep = ""))
    if (sum(is.na(data[, id])) > 0) stop(paste("Missing values in variable \"", id, "\".", sep = ""))
    if (sum(is.na(data[, time])) > 0) stop(paste("Missing values in variable \"", time, "\".", sep = ""))

    data <- data[order(data[, id], data[, time]), ]
    Ttot <- sort(unique(data[, time]))
    units <- unique(data[, id])

    if (nrow(data) != TT * N) {
        stop("Data are not balanced or \"index\" does not uniquely identity an observation.")
    }

    D.sav <- D <- matrix(data[, Dname], TT, N)
    D <- apply(D, 2, function(vec) cumsum(vec))
    T0 <- TT - D[TT, ]

    id.drop <- which(T0 <= trim.npre)
    N.drop <- length(id.drop)
    D <- ifelse(D > 0, 1, 0)
    if (sum(abs(D - D.sav)) != 0) {
        cat("\nTreatment status changed to \"treated\" after a unit has ever been treated; no switch on-and-off allowed.\n")
    }
    if (N.drop > 0) {
        N <- N - N.drop
        D <- D[, -id.drop, drop = FALSE]
        data <- data[rep(T0, each = TT) > trim.npre, ]
        units <- units[-id.drop]
        T0 <- T0[-id.drop]
        cat(paste0("\nDrop ", length(id.drop), " units with ", trim.npre, " or fewer pre-treatment periods.\n"))
    }

    treat <- ifelse(D[TT, ] == 1, 1, 0)
    id.tr <- which(treat == 1)
    id.co <- which(treat == 0)
    Ntr <- length(id.tr)
    Nco <- length(id.co)
    if (Ntr == 0) stop("No treated units remain.")
    if (Nco == 0) stop("No control units remain.")

    if (Ntr <= 5) {
        cat("Too few treated unit(s). Uncertainty estimates not provided.\n")
        vce <- "none"
    }

    T0.tr <- T0[id.tr]
    T0.min <- min(T0.tr)

    if (Ntr == 1) {
        sameT0 <- TRUE
    } else {
        sameT0 <- (var(T0.tr) == 0)
    }
    if (sameT0) Tpre <- Ttot[1:unique(T0.tr)]

    outcome <- matrix(data[, Yname], N, TT, byrow = TRUE)
    Y.var <- paste0(Yname, Ttot)
    colnames(outcome) <- Y.var

    if (!is.factor(data[, id])) data[, id] <- as.factor(data[, id])

    if (p > 0) {
        if (!is.null(X.avg.time)) {
            if (sameT0 == FALSE) {
                stop("\"X.avg.time\" is only allowed when the treatment starts at the same time.")
            }
            if (is.list(X.avg.time)) {
                if (length(X.avg.time) != p) {
                    stop("Length of \"X.avg.time\" (as a list) must equal the number of covariates.")
                }
                Xvar <- matrix(NA, N, p)
                colnames(Xvar) <- Xname
                for (i in 1:p) {
                    this.period <- X.avg.time[[i]]
                    if (sum(1 - this.period %in% Tpre) > 0) {
                        stop("Elements in \"X.avg.time\" must be in the pre-treatment period.")
                    }
                    selected.row <- which(data[, time] %in% this.period)
                    X.pre <- data[selected.row, c(id, Xname[i]), drop = FALSE]
                    agg <- aggregate(X.pre[, Xname[i], drop = FALSE],
                                     list(unit = X.pre[, id]), mean, na.rm = TRUE)
                    covar.tmp <- agg[, -1]
                    if (length(covar.tmp) != N) {
                        stop(paste0("Missing values in ", Xname[i], " in specified years."))
                    } else {
                        Xvar[, i] <- covar.tmp
                    }
                }
            } else {
                if (sum(1 - X.avg.time %in% Tpre) > 0) {
                    stop("\"X.avg.time\" must be in the pre-treatment period.")
                }
                selected.row <- which(data[, time] %in% X.avg.time)
                X.pre <- data[selected.row, Xname, drop = FALSE]
                agg <- aggregate(X.pre, list(unit = data[selected.row, id]), mean, na.rm = TRUE)
                Xvar <- agg[, -1, drop = FALSE]
                if (nrow(Xvar) != N) stop("Missing values in covariates.")
            }
            for (i in 1:p) {
                if (sum(is.na(Xvar[, i])) > 0) {
                    stop(paste0("Missing values in variable \"", Xname[i], "\".", sep = ""))
                }
            }
        } else {
            Xvar <- matrix(NA, N, p)
            colnames(Xvar) <- Xname
            for (i in 1:p) {
                if (sum(is.na(data[, Xname[i]])) > 0) {
                    warning(paste0("Missing values in variable \"", Xname[i], "\".", sep = ""))
                }
                X.tmp <- matrix(data[, Xname[i]], N, TT, byrow = TRUE)
                X.var <- apply(X.tmp, 1, var, na.rm = TRUE)
                if (sum(is.na(X.var)) > 0) {
                    stop(paste0("Variable \"", Xname[i], "\" is completely missing in some unit(s)."))
                }
                if (sum(X.var) != 0) {
                    stop(paste0("\"", Xname[i], "\" is not time-invariant for some unit(s)."))
                }
                Xvar[, i] <- apply(X.tmp, 1, mean, na.rm = TRUE)
            }
        }
    }

    if (p > 0) {
        data.wide <- cbind.data.frame(id = 1:N, unit = units, treat = treat, T0 = T0, outcome, Xvar)
    } else {
        data.wide <- cbind.data.frame(id = 1:N, unit = units, treat = treat, T0 = T0, outcome)
    }

    if (sameT0) {
        bal.out <- tjbal_estimate_single(
            data = data.wide, Y = Yname, D = "treat", X = Xname,
            Y.match.time = Y.match.time, Y.match.npre = Y.match.npre,
            Ttot = Ttot,
            demean = demean, estimator = estimator, sigma = sigma,
            eta = eta, nystrom = nystrom, nystrom_m = nystrom_m,
            print.baltable = print.baltable,
            vce = vce, conf.lvl = conf.lvl,
            nsims = nsims, parallel = parallel, cores = cores
        )
    } else {
        bal.out <- tjbal_estimate_multi(
            data = data.wide, Y = Yname, D = "treat", X = Xname,
            Y.match.time = Y.match.time, Y.match.npre = Y.match.npre,
            Ttot = Ttot,
            demean = demean, estimator = estimator, sigma = sigma,
            eta = eta, nystrom = nystrom, nystrom_m = nystrom_m,
            vce = vce, conf.lvl = conf.lvl,
            nsims = nsims, parallel = parallel, cores = cores
        )
    }

    out <- c(list(sameT0 = sameT0, index = index, Yname = Yname), bal.out)
    out$call <- match.call()
    class(out) <- "tjbal"
    return(out)
}
