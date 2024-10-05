> DESeq2:::fitNbinomGLMs
function (object, modelMatrix = NULL, modelFormula, alpha_hat,
    lambda, renameCols = TRUE, betaTol = 1e-08, maxit = 100,
    useOptim = TRUE, useQR = TRUE, forceOptim = FALSE, warnNonposVar = TRUE,
    minmu = 0.5, type = c("DESeq2", "glmGamPoi"))
{
    type <- match.arg(type, c("DESeq2", "glmGamPoi"))
    if (missing(modelFormula)) {
        modelFormula <- design(object)
    }
    if (is.null(modelMatrix)) {
        modelAsFormula <- TRUE
        modelMatrix <- stats::model.matrix.default(modelFormula,
            data = as.data.frame(colData(object)))
    }
    else {
        modelAsFormula <- FALSE
    }
    stopifnot(all(MatrixGenerics::colSums(abs(modelMatrix)) >
        0))
    modelMatrixNames <- colnames(modelMatrix)
    modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
    modelMatrixNames <- make.names(modelMatrixNames)
    if (renameCols) {
        convertNames <- DESeq2:::renameModelMatrixColumns(colData(object),
            modelFormula)
        convertNames <- convertNames[convertNames$from %in% modelMatrixNames,
            , drop = FALSE]
        modelMatrixNames[match(convertNames$from, modelMatrixNames)] <- convertNames$to
    }
    colnames(modelMatrix) <- modelMatrixNames
    normalizationFactors <- DESeq2:::getSizeOrNormFactors(object)
    if (missing(alpha_hat)) {
        alpha_hat <- dispersions(object)
    }
    if (length(alpha_hat) != nrow(object)) {
        stop("alpha_hat needs to be the same length as nrows(object)")
    }
    if (missing(lambda)) {
        lambda <- rep(1e-06, ncol(modelMatrix))
    }
    wlist <- DESeq2:::getAndCheckWeights(object, modelMatrix)
    weights <- wlist$weights
    useWeights <- wlist$useWeights
    if (type == "glmGamPoi") {
        stopifnot(`type = 'glmGamPoi' cannot handle weights` = !useWeights,
            `type = 'glmGamPoi' does not support NA's in alpha_hat` = all(!is.na(alpha_hat)))
        gp_res <- glmGamPoi::glm_gp(counts(object), design = modelMatrix,
            size_factors = FALSE, offset = log(normalizationFactors),
            overdispersion = alpha_hat, verbose = FALSE)
        logLikeMat <- dnbinom(counts(object), mu = gp_res$Mu,
            size = 1/alpha_hat, log = TRUE)
        logLike <- MatrixGenerics::rowSums(logLikeMat)
        res <- list(logLike = logLike, betaConv = rep(TRUE, nrow(object)),
            betaMatrix = gp_res$Beta/log(2), betaSE = NULL, mu = gp_res$Mu,
            betaIter = rep(NA, nrow(object)), modelMatrix = modelMatrix,
            nterms = ncol(modelMatrix), hat_diagonals = NULL)
        return(res)
    }
    justIntercept <- if (modelAsFormula) {
        modelFormula == formula(~1)
    }else {
        ncol(modelMatrix) == 1 & all(modelMatrix == 1)
    }
    if (justIntercept & all(lambda <= 1e-06)) {
        alpha <- alpha_hat
        betaConv <- rep(TRUE, nrow(object))
        betaIter <- rep(1, nrow(object))
        betaMatrix <- if (useWeights) {
            matrix(log2(MatrixGenerics::rowSums(weights * counts(object,
                normalized = TRUE))/MatrixGenerics::rowSums(weights)),
                ncol = 1)
        }
        else {
            matrix(log2(MatrixGenerics::rowMeans(counts(object,
                normalized = TRUE))), ncol = 1)
        }
        mu <- normalizationFactors * as.numeric(2^betaMatrix)
        logLikeMat <- dnbinom(counts(object), mu = mu, size = 1/alpha,
            log = TRUE)
        logLike <- if (useWeights) {
            MatrixGenerics::rowSums(weights * logLikeMat)
        } else {
            MatrixGenerics::rowSums(logLikeMat)
        }
        modelMatrix <- stats::model.matrix.default(~1, data = as.data.frame(colData(object)))
        colnames(modelMatrix) <- modelMatrixNames <- "Intercept"
        w <- if (useWeights) {
            weights * (mu^-1 + alpha)^-1
        }else {
            (mu^-1 + alpha)^-1
        }
        xtwx <- MatrixGenerics::rowSums(w)
        sigma <- xtwx^-1
        betaSE <- matrix(log2(exp(1)) * sqrt(sigma), ncol = 1)
        hat_diagonals <- w * xtwx^-1
        res <- list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
            betaSE = betaSE, mu = mu, betaIter = betaIter, modelMatrix = modelMatrix,
            nterms = 1, hat_diagonals = hat_diagonals)
        return(res)
    }
    qrx <- qr(modelMatrix)
    if (qrx$rank == ncol(modelMatrix)) {
        Q <- qr.Q(qrx)
        R <- qr.R(qrx)
        y <- t(log(counts(object, normalized = TRUE) + 0.1))
        beta_mat <- t(solve(R, t(Q) %*% y))
    }
    else {
        if ("Intercept" %in% modelMatrixNames) {
            beta_mat <- matrix(0, ncol = ncol(modelMatrix), nrow = nrow(object))
            logBaseMean <- log(MatrixGenerics::rowMeans(counts(object,
                normalized = TRUE)))
            beta_mat[, which(modelMatrixNames == "Intercept")] <- logBaseMean
        }
        else {
            beta_mat <- matrix(1, ncol = ncol(modelMatrix), nrow = nrow(object))
        }
    }
    lambdaNatLogScale <- lambda/log(2)^2
    betaRes <- DESeq2:::fitBetaWrapper(ySEXP = counts(object), xSEXP = modelMatrix,
        nfSEXP = normalizationFactors, alpha_hatSEXP = alpha_hat,
        beta_matSEXP = beta_mat, lambdaSEXP = lambdaNatLogScale,
        weightsSEXP = weights, useWeightsSEXP = useWeights, tolSEXP = betaTol,
        maxitSEXP = maxit, useQRSEXP = useQR, minmuSEXP = minmu)
    mu <- normalizationFactors * t(exp(modelMatrix %*% t(betaRes$beta_mat)))
    dispersionVector <- rep(dispersions(object), times = ncol(object))
    logLike <- nbinomLogLike(counts(object), mu, dispersions(object),
        weights, useWeights)
    rowStable <- apply(betaRes$beta_mat, 1, function(row) sum(is.na(row))) ==
        0
    rowVarPositive <- apply(betaRes$beta_var_mat, 1, function(row) sum(row <=
        0)) == 0
    betaConv <- betaRes$iter < maxit
    betaMatrix <- log2(exp(1)) * betaRes$beta_mat
    colnames(betaMatrix) <- modelMatrixNames
    colnames(modelMatrix) <- modelMatrixNames
    betaSE <- log2(exp(1)) * sqrt(pmax(betaRes$beta_var_mat,
        0))
    colnames(betaSE) <- paste0("SE_", modelMatrixNames)
    rowsForOptim <- if (useOptim) {
        which(!betaConv | !rowStable | !rowVarPositive)
    }
    else {
        which(!rowStable | !rowVarPositive)
    }
    if (forceOptim) {
        rowsForOptim <- seq_along(betaConv)
    }
    if (length(rowsForOptim) > 0) {
        resOptim <- fitNbinomGLMsOptim(object, modelMatrix, lambda,
            rowsForOptim, rowStable, normalizationFactors, alpha_hat,
            weights, useWeights, betaMatrix, betaSE, betaConv,
            beta_mat, mu, logLike, minmu = minmu)
        betaMatrix <- resOptim$betaMatrix
        betaSE <- resOptim$betaSE
        betaConv <- resOptim$betaConv
        mu <- resOptim$mu
        logLike <- resOptim$logLike
    }
    stopifnot(!any(is.na(betaSE)))
    nNonposVar <- sum(MatrixGenerics::rowSums(betaSE == 0) >
        0)
    if (warnNonposVar & nNonposVar > 0)
        warning(nNonposVar, "rows had non-positive estimates of variance for coefficients")
    list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
        betaSE = betaSE, mu = mu, betaIter = betaRes$iter, modelMatrix = modelMatrix,
        nterms = ncol(modelMatrix), hat_diagonals = betaRes$hat_diagonals)
}


> DESeq2:::fitNbinomGLMsOptim
function (object, modelMatrix, lambda, rowsForOptim, rowStable,
    normalizationFactors, alpha_hat, weights, useWeights, betaMatrix,
    betaSE, betaConv, beta_mat, mu, logLike, minmu = 0.5)
{
    x <- modelMatrix
    lambdaNatLogScale <- lambda/log(2)^2
    large <- 30
    for (row in rowsForOptim) {
        betaRow <- if (rowStable[row] & all(abs(betaMatrix[row,
            ]) < large)) {
            betaMatrix[row, ]
        }
        else {
            beta_mat[row, ]
        }
        nf <- normalizationFactors[row, ]
        k <- counts(object)[row, ]
        alpha <- alpha_hat[row]
        objectiveFn <- function(p) {
            mu_row <- as.numeric(nf * 2^(x %*% p))
            logLikeVector <- dnbinom(k, mu = mu_row, size = 1/alpha,
                log = TRUE)
            logLike <- if (useWeights) {
                sum(weights[row, ] * logLikeVector)
            }
            else {
                sum(logLikeVector)
            }
            logPrior <- sum(dnorm(p, 0, sqrt(1/lambda), log = TRUE))
            negLogPost <- -1 * (logLike + logPrior)
            if (is.finite(negLogPost))
                negLogPost
            else 10^300
        }
        o <- optim(betaRow, objectiveFn, method = "L-BFGS-B",
            lower = -large, upper = large)
        ridge <- if (length(lambdaNatLogScale) > 1) {
            diag(lambdaNatLogScale)
        }
        else {
            as.matrix(lambdaNatLogScale, ncol = 1)
        }
        if (o$convergence == 0) {
            betaConv[row] <- TRUE
        }
        betaMatrix[row, ] <- o$par
        mu_row <- as.numeric(nf * 2^(x %*% o$par))
        mu[row, ] <- mu_row
        mu_row[mu_row < minmu] <- minmu
        w <- if (useWeights) {
            diag(weights[row, ] * (mu_row^-1 + alpha)^-1)
        }
        else {
            diag((mu_row^-1 + alpha)^-1)
        }
        xtwx <- t(x) %*% w %*% x
        xtwxRidgeInv <- solve(xtwx + ridge)
        sigma <- xtwxRidgeInv %*% xtwx %*% xtwxRidgeInv
        betaSE[row, ] <- log2(exp(1)) * sqrt(pmax(diag(sigma),
            0))
        logLikeVector <- dnbinom(k, mu = mu_row, size = 1/alpha,
            log = TRUE)
        logLike[row] <- if (useWeights) {
            sum(weights[row, ] * logLikeVector)
        }
        else {
            sum(logLikeVector)
        }
    }
    return(list(betaMatrix = betaMatrix, betaSE = betaSE, betaConv = betaConv,
        mu = mu, logLike = logLike))
}

 glm.fitter
function (x, y, weights = rep.int(1, nobs), start = NULL, etastart = NULL,
    mustart = NULL, offset = rep.int(0, nobs), family = gaussian(),
    control = list(), intercept = TRUE, singular.ok = TRUE)
{
    control <- do.call("glm.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    valideta <- family$valideta %||% function(eta) TRUE
    validmu <- family$validmu %||% function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta))
            stop("invalid linear predictor values in empty model",
                call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu))
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- sqrt((weights * mu.eta(eta)^2)/variance(mu))
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep_len(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    }
    else {
        coefold <- NULL
        eta <- etastart %||% {
            if (!is.null(start))
                if (length(start) != nvars)
                  stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                    nvars, paste(deparse(xnames), collapse = ", ")),
                    domain = NA)
                else {
                  coefold <- start
                  offset + as.vector(if (NCOL(x) == 1L)
                    x * start
                  else x %*% start)
                }
            else family$linkfun(mustart)
        }
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some",
                call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
        for (iter in 1L:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (anyNA(varmu))
                stop("NAs in V(mu)")
            if (any(varmu == 0))
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good])))
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning(gettextf("no observations informative at iteration %d",
                  iter), domain = NA)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
            fit <- .Call(C_Cdqrls, x[good, , drop = FALSE] *
                w, z * w, min(1e-07, control$epsilon/1000), check = FALSE)
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d",
                  iter), domain = NA)
                break
            }
            if (nobs < fit$rank)
                stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation",
                  "X matrix has rank %d, but only %d observations"),
                  fit$rank, nobs), domain = NA)
            if (!singular.ok && fit$rank < nvars)
                stop("singular fit encountered")
            start[fit$pivot] <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace)
                cat("Deviance = ", dev, " Iterations - ", iter,
                  "\n", sep = "")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                warning("step size truncated due to divergence",
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit)
                    stop("inner loop 1; cannot correct step size",
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace)
                  cat("Step halved: new deviance = ", dev, "\n",
                    sep = "")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                warning("step size truncated: out of bounds",
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit)
                    stop("inner loop 2; cannot correct step size",
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace)
                  cat("Step halved: new deviance = ", dev, "\n",
                    sep = "")
            }
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {
                devold <- dev
                coef <- coefold <- start
            }
        }
        if (!conv)
            warning("glm.fit: algorithm did not converge", call. = FALSE)
        if (boundary)
            warning("glm.fit: algorithm stopped at boundary value",
                call. = FALSE)
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps))
                warning("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                  call. = FALSE)
        }
        if (family$family == "poisson") {
            if (any(mu < eps))
                warning("glm.fit: fitted rates numerically 0 occurred",
                  call. = FALSE)
        }
        if (fit$rank < nvars)
            coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        xxnames <- xnames[fit$pivot]
        residuals <- (y - mu)/mu.eta(eta)
        fit$qr <- as.matrix(fit$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
        }
        else Rmat <- fit$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY)
        names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("",
            sum(good) - fit$rank))
    wtdmu <- if (intercept)
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY)
        0
    else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu,
        effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat,
        rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank",
            "qraux", "pivot", "tol")], class = "qr"), family = family,
        linear.predictors = eta, deviance = dev, aic = aic.model,
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,
        df.residual = resdf, df.null = nulldf, y = y, converged = conv,
        boundary = boundary)
}


> glm_gp
function (data, design = ~1, col_data = NULL, reference_level = NULL,
    offset = 0, size_factors = c("normed_sum", "deconvolution",
        "poscounts", "ratio"), overdispersion = TRUE, overdispersion_shrinkage = TRUE,
    ridge_penalty = 0, do_cox_reid_adjustment = TRUE, subsample = FALSE,
    on_disk = NULL, use_assay = NULL, verbose = FALSE)
{
    if (inherits(data, "formula")) {
        if (length(design) != 2 || design != ~1) {
            stop("If the first argument is already a formula, the second argument must not be set. Please call this function like this:\n",
                "'glm_gp(data = mat, design = ~ a + b + c, ...)'",
                call. = FALSE)
        }
        extr <- extract_data_from_formula(data, col_data, parent.frame())
        data <- extr$data
        design <- extr$design
    }
    if (is.vector(data)) {
        data <- matrix(data, nrow = 1)
    }
    data_mat <- handle_data_parameter(data, use_assay, on_disk,
        verbose = verbose)
    col_data <- get_col_data(data, col_data)
    des <- handle_design_parameter(design, data, col_data, reference_level)
    res <- glm_gp_impl(data_mat, model_matrix = des$model_matrix,
        offset = offset, size_factors = size_factors, overdispersion = overdispersion,
        overdispersion_shrinkage = overdispersion_shrinkage,
        ridge_penalty = ridge_penalty, do_cox_reid_adjustment = do_cox_reid_adjustment,
        subsample = subsample, verbose = verbose)
    rownames(data_mat) <- rownames(data)
    colnames(data_mat) <- colnames(data)
    res$data <- SummarizedExperiment::SummarizedExperiment(list(counts = data_mat),
        colData = col_data)
    res$model_matrix <- des$model_matrix
    if (is.null(colnames(res$model_matrix))) {
        colnames(res$model_matrix) <- paste0("Coef_", seq_len(ncol(res$model_matrix)))
    }
    res$design_formula <- des$design_formula
    colnames(res$Beta) <- colnames(res$model_matrix)
    rownames(res$Beta) <- rownames(data)
    if (!is.null(res$ridge_penalty)) {
        names(res$ridge_penalty) <- colnames(res$model_matrix)
    }
    rownames(res$Mu) <- rownames(data)
    colnames(res$Mu) <- colnames(data)
    rownames(res$Offset) <- rownames(data)
    colnames(res$Offset) <- colnames(data)
    names(res$overdispersions) <- rownames(data)
    names(res$deviances) <- rownames(data)
    names(res$size_factors) <- colnames(data)
    class(res) <- "glmGamPoi"
    res
}



> DESeq2:::fitBetaWrapper
function (ySEXP, xSEXP, nfSEXP, alpha_hatSEXP, contrastSEXP,
    beta_matSEXP, lambdaSEXP, weightsSEXP, useWeightsSEXP, tolSEXP,
    maxitSEXP, useQRSEXP, minmuSEXP)
{
    if (missing(contrastSEXP)) {
        contrastSEXP <- c(1, rep(0, ncol(xSEXP) - 1))
    }
    arg.names <- names(formals(fitBetaWrapper))
    na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
    if (any(na.test))
        stop(paste("in call to fitBeta, the following arguments contain NA:",
            paste(arg.names[na.test], collapse = ", ")))
    fitBeta(ySEXP = ySEXP, xSEXP = xSEXP, nfSEXP = nfSEXP, alpha_hatSEXP = alpha_hatSEXP,
        contrastSEXP = contrastSEXP, beta_matSEXP = beta_matSEXP,
        lambdaSEXP = lambdaSEXP, weightsSEXP = weightsSEXP, useWeightsSEXP = useWeightsSEXP,
        tolSEXP = tolSEXP, maxitSEXP = maxitSEXP, useQRSEXP = useQRSEXP,
        minmuSEXP = minmuSEXP)
}

> DESeq2:::fitNbinomGLMs # This is the function that fits the negative binomial GLMs, estimages the coefficients and the standard errors
function (object, modelMatrix = NULL, modelFormula, alpha_hat,
    lambda, renameCols = TRUE, betaTol = 1e-08, maxit = 100,
    useOptim = TRUE, useQR = TRUE, forceOptim = FALSE, warnNonposVar = TRUE,
    minmu = 0.5, type = c("DESeq2", "glmGamPoi"))
{
    type <- match.arg(type, c("DESeq2", "glmGamPoi"))
    if (missing(modelFormula)) {
        modelFormula <- design(object)
    }
    if (is.null(modelMatrix)) {
        modelAsFormula <- TRUE
        modelMatrix <- stats::model.matrix.default(modelFormula,
            data = as.data.frame(colData(object)))
    } else {
        modelAsFormula <- FALSE
    }
    stopifnot(all(MatrixGenerics::colSums(abs(modelMatrix)) >
        0))
    modelMatrixNames <- colnames(modelMatrix)
    modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
    modelMatrixNames <- make.names(modelMatrixNames)
    if (renameCols) {
        convertNames <- DESeq2:::renameModelMatrixColumns(colData(object),
            modelFormula)
        convertNames <- convertNames[convertNames$from %in% modelMatrixNames,
            , drop = FALSE]
        modelMatrixNames[match(convertNames$from, modelMatrixNames)] <- convertNames$to
    }
    colnames(modelMatrix) <- modelMatrixNames
    normalizationFactors <- DESeq2:::getSizeOrNormFactors(object)
    if (missing(alpha_hat)) {
        alpha_hat <- dispersions(object)
    }
    if (length(alpha_hat) != nrow(object)) {
        stop("alpha_hat needs to be the same length as nrows(object)")
    }
    if (missing(lambda)) {
        lambda <- rep(1e-06, ncol(modelMatrix))
    }
    wlist <- DESeq2:::getAndCheckWeights(object, modelMatrix)
    weights <- wlist$weights
    useWeights <- wlist$useWeights
    if (type == "glmGamPoi") {
        stopifnot(`type = 'glmGamPoi' cannot handle weights` = !useWeights,
            `type = 'glmGamPoi' does not support NA's in alpha_hat` = all(!is.na(alpha_hat)))
        gp_res <- glmGamPoi::glm_gp(counts(object), design = modelMatrix,
            size_factors = FALSE, offset = log(normalizationFactors),
            overdispersion = alpha_hat, verbose = FALSE)
        logLikeMat <- dnbinom(counts(object), mu = gp_res$Mu,
            size = 1/alpha_hat, log = TRUE)
        logLike <- MatrixGenerics::rowSums(logLikeMat)
        res <- list(logLike = logLike, betaConv = rep(TRUE, nrow(object)),
            betaMatrix = gp_res$Beta/log(2), betaSE = NULL, mu = gp_res$Mu,
            betaIter = rep(NA, nrow(object)), modelMatrix = modelMatrix,
            nterms = ncol(modelMatrix), hat_diagonals = NULL)
        return(res)
    }
    justIntercept <- if (modelAsFormula) {
        modelFormula == formula(~1)
    } else {
        ncol(modelMatrix) == 1 & all(modelMatrix == 1)
    }
    if (justIntercept & all(lambda <= 1e-06)) {
        alpha <- alpha_hat
        betaConv <- rep(TRUE, nrow(object))
        betaIter <- rep(1, nrow(object))
        betaMatrix <- if (useWeights) {
            matrix(log2(MatrixGenerics::rowSums(weights * counts(object,
                normalized = TRUE))/MatrixGenerics::rowSums(weights)),
                ncol = 1)
        }else {    matrix(log2(MatrixGenerics::rowMeans(counts(object,
                normalized = TRUE))), ncol = 1)
        }
        mu <- normalizationFactors * as.numeric(2^betaMatrix)
        logLikeMat <- dnbinom(counts(object), mu = mu, size = 1/alpha,
            log = TRUE)
        logLike <- if (useWeights) {
            MatrixGenerics::rowSums(weights * logLikeMat)
        }else {
            MatrixGenerics::rowSums(logLikeMat)
        }
        modelMatrix <- stats::model.matrix.default(~1, data = as.data.frame(colData(object)))
        colnames(modelMatrix) <- modelMatrixNames <- "Intercept"
        w <- if (useWeights) {
            weights * (mu^-1 + alpha)^-1
        }else {
            (mu^-1 + alpha)^-1
        }
        xtwx <- MatrixGenerics::rowSums(w)
        sigma <- xtwx^-1
        betaSE <- matrix(log2(exp(1)) * sqrt(sigma), ncol = 1)
        hat_diagonals <- w * xtwx^-1
        res <- list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
            betaSE = betaSE, mu = mu, betaIter = betaIter, modelMatrix = modelMatrix,
            nterms = 1, hat_diagonals = hat_diagonals)
        return(res)
    }
    qrx <- qr(modelMatrix)
    if (qrx$rank == ncol(modelMatrix)) {
        Q <- qr.Q(qrx)
        R <- qr.R(qrx)
        y <- t(log(counts(object, normalized = TRUE) + 0.1))
        beta_mat <- t(solve(R, t(Q) %*% y))
    } else {
        if ("Intercept" %in% modelMatrixNames) {
            beta_mat <- matrix(0, ncol = ncol(modelMatrix), nrow = nrow(object))
            logBaseMean <- log(MatrixGenerics::rowMeans(counts(object,
                normalized = TRUE)))
            beta_mat[, which(modelMatrixNames == "Intercept")] <- logBaseMean
        }   else {
            beta_mat <- matrix(1, ncol = ncol(modelMatrix), nrow = nrow(object))
        }
    }
    lambdaNatLogScale <- lambda/log(2)^2
    betaRes <- fitBetaWrapper(ySEXP = counts(object), xSEXP = modelMatrix,
        nfSEXP = normalizationFactors, alpha_hatSEXP = alpha_hat,
        beta_matSEXP = beta_mat, lambdaSEXP = lambdaNatLogScale,
        weightsSEXP = weights, useWeightsSEXP = useWeights, tolSEXP = betaTol,
        maxitSEXP = maxit, useQRSEXP = useQR, minmuSEXP = minmu)
    mu <- normalizationFactors * t(exp(modelMatrix %*% t(betaRes$beta_mat)))
    dispersionVector <- rep(dispersions(object), times = ncol(object))
    logLike <- DESeq2:::nbinomLogLike(counts(object), mu, dispersions(object),
        weights, useWeights)
    rowStable <- apply(betaRes$beta_mat, 1, function(row) sum(is.na(row))) ==
        0
    rowVarPositive <- apply(betaRes$beta_var_mat, 1, function(row) sum(row <=
        0)) == 0
    betaConv <- betaRes$iter < maxit
    betaMatrix <- log2(exp(1)) * betaRes$beta_mat # convert to log2 fold change
    colnames(betaMatrix) <- modelMatrixNames
    colnames(modelMatrix) <- modelMatrixNames
    betaSE <- log2(exp(1)) * sqrt(pmax(betaRes$beta_var_mat,
        0))
    colnames(betaSE) <- paste0("SE_", modelMatrixNames)
    rowsForOptim <- if (useOptim) {
        which(!betaConv | !rowStable | !rowVarPositive)
    }else {
        which(!rowStable | !rowVarPositive)
    }
    if (forceOptim) {
        rowsForOptim <- seq_along(betaConv)
    }
    if (length(rowsForOptim) > 0) {
        resOptim <- DESeq2:::fitNbinomGLMsOptim(object, modelMatrix, lambda,
            rowsForOptim, rowStable, normalizationFactors, alpha_hat,
            weights, useWeights, betaMatrix, betaSE, betaConv,
            beta_mat, mu, logLike, minmu = minmu)
        betaMatrix <- resOptim$betaMatrix
        betaSE <- resOptim$betaSE
        betaConv <- resOptim$betaConv
        mu <- resOptim$mu
        logLike <- resOptim$logLike
    }
    stopifnot(!any(is.na(betaSE)))
    nNonposVar <- sum(MatrixGenerics::rowSums(betaSE == 0) >
        0)
    if (warnNonposVar & nNonposVar > 0)
        warning(nNonposVar, "rows had non-positive estimates of variance for coefficients")
    list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
        betaSE = betaSE, mu = mu, betaIter = betaRes$iter, modelMatrix = modelMatrix,
        nterms = ncol(modelMatrix), hat_diagonals = betaRes$hat_diagonals)
}
<bytecode: 0x7f934d8fcba

# Load necessary libraries
library(DESeq2)

# Simulated count data and metadata
counts <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
colData <- data.frame(condition = factor(rep(c("control", "treatment"), each = 5)))

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)

# Run the DESeq2 analysis pipeline
dds <- DESeq(dds)

# 1. Extract the design matrix
modelMatrix <- model.matrix(design(dds), colData(dds))

# 2. Extract the beta coefficients (log2 fold changes)
betaMatrix <- coef(dds)

# 3. Extract the size factors
sizeFactors <- sizeFactors(dds)

# 4. Calculate the fitted mean counts (mu_ij)
fitted_values <- exp(modelMatrix %*% t(betaMatrix)) * sizeFactors

# View the fitted values (mean counts)
head(fitted_values)



> nbinomWaldTest
function (object, betaPrior = FALSE, betaPriorVar, modelMatrix = NULL,
    modelMatrixType, betaTol = 1e-08, maxit = 100, useOptim = TRUE,
    quiet = FALSE, useT = FALSE, df, useQR = TRUE, minmu = 0.5)
{
    if (is.null(dispersions(object))) {
        stop("testing requires dispersion estimates, first call estimateDispersions()")
    }
    stopifnot(length(maxit) == 1)
    object <- DESeq2:::sanitizeRowRanges(object)
    if ("results" %in% mcols(mcols(object))$type) {
        if (!quiet)
            message("found results columns, replacing these")
        object <- removeResults(object)
    }
    if (is.null(mcols(object)$allZero)) {
        object <- getBaseMeansAndVariances(object)
    }
    objectNZ <- object[!mcols(object)$allZero, , drop = FALSE]
    if (is.null(modelMatrix)) {
        modelAsFormula <- TRUE
        termsOrder <- attr(terms.formula(design(object)), "order")
        interactionPresent <- any(termsOrder > 1)
        if (missing(betaPrior)) {
            betaPrior <- FALSE
        }
        DESeq2:::designAndArgChecker(object, betaPrior)
        stopifnot(is.logical(betaPrior))
        blindDesign <- design(object) == formula(~1)
        if (blindDesign) {
            betaPrior <- FALSE
        }
        if (missing(modelMatrixType) || is.null(modelMatrixType)) {
            modelMatrixType <- if (betaPrior) {
                "expanded"
            }
            else {
                "standard"
            }
        }
        if (modelMatrixType == "expanded" & !betaPrior) {
            stop("expanded model matrices require a beta prior")
        }
        attr(object, "modelMatrixType") <- modelMatrixType
        hasIntercept <- attr(terms(design(object)), "intercept") ==
            1
        renameCols <- hasIntercept
    }
    else {
        if (missing(betaPrior)) {
            betaPrior <- FALSE
        }
        if (betaPrior) {
            if (missing(betaPriorVar))
                stop("user-supplied model matrix with betaPrior=TRUE requires supplying betaPriorVar")
        }
        modelAsFormula <- FALSE
        attr(object, "modelMatrixType") <- "user-supplied"
        renameCols <- FALSE
    }
    if (!betaPrior) {
        fit <- fitNbinomGLMs(objectNZ, betaTol = betaTol, maxit = maxit,
            useOptim = useOptim, useQR = useQR, renameCols = renameCols,
            modelMatrix = modelMatrix, minmu = minmu)
        H <- fit$hat_diagonals
        mu <- fit$mu
        modelMatrix <- fit$modelMatrix
        modelMatrixNames <- fit$modelMatrixNames
        betaPriorVar <- rep(1e+06, ncol(fit$modelMatrix))
    }
    else {
        priorFitList <- fitGLMsWithPrior(object = object, betaTol = betaTol,
            maxit = maxit, useOptim = useOptim, useQR = useQR,
            betaPriorVar = betaPriorVar, modelMatrix = modelMatrix,
            minmu = minmu)
        fit <- priorFitList$fit
        H <- priorFitList$H
        mu <- priorFitList$mu
        betaPriorVar <- priorFitList$betaPriorVar
        modelMatrix <- priorFitList$modelMatrix
        mleBetaMatrix <- priorFitList$mleBetaMatrix
        mcols(object) <- mcols(object)[, grep("MLE_", names(mcols(object)),
            invert = TRUE)]
    }
    dimnames(mu) <- NULL
    assays(objectNZ, withDimnames = FALSE)[["mu"]] <- mu
    assays(object, withDimnames = FALSE)[["mu"]] <- buildMatrixWithNARows(mu,
        mcols(object)$allZero)
    dimnames(H) <- NULL
    assays(objectNZ, withDimnames = FALSE)[["H"]] <- H
    assays(object, withDimnames = FALSE)[["H"]] <- buildMatrixWithNARows(H,
        mcols(object)$allZero)
    attr(object, "betaPrior") <- betaPrior
    attr(object, "betaPriorVar") <- betaPriorVar
    attr(object, "modelMatrix") <- modelMatrix
    attr(object, "test") <- "Wald"
    dispModelMatrix <- if (modelAsFormula) {
        getModelMatrix(object)
    }
    else {
        modelMatrix
    }
    attr(object, "dispModelMatrix") <- dispModelMatrix
    cooks <- calculateCooksDistance(objectNZ, H, dispModelMatrix)
    maxCooks <- recordMaxCooks(design(object), colData(object),
        dispModelMatrix, cooks, nrow(objectNZ))
    assays(object, withDimnames = FALSE)[["cooks"]] <- buildMatrixWithNARows(cooks,
        mcols(object)$allZero)
    modelMatrixNames <- colnames(modelMatrix)
    betaMatrix <- fit$betaMatrix
    colnames(betaMatrix) <- modelMatrixNames
    betaSE <- fit$betaSE
    colnames(betaSE) <- paste0("SE_", modelMatrixNames)
    WaldStatistic <- betaMatrix/betaSE
    colnames(WaldStatistic) <- paste0("WaldStatistic_", modelMatrixNames)
    if (useT) {
        if (!missing(df)) {
            stopifnot(length(df) == 1 | length(df) == nrow(object))
            if (length(df) == 1) {
                df <- rep(df, nrow(objectNZ))
            }
            else {
                df <- df[!mcols(object)$allZero]
            }
        }
        else {
            if ("weights" %in% assayNames(object)) {
                wlist <- getAndCheckWeights(objectNZ, dispModelMatrix)
                num.samps <- MatrixGenerics::rowSums(wlist$weights)
            }
            else {
                num.samps <- rep(ncol(object), nrow(objectNZ))
            }
            df <- num.samps - ncol(dispModelMatrix)
        }
        df <- ifelse(df > 0, df, NA)
        stopifnot(length(df) == nrow(WaldStatistic))
        WaldPvalue <- 2 * pt(abs(WaldStatistic), df = df, lower.tail = FALSE)
    }
    else {
        WaldPvalue <- 2 * pnorm(abs(WaldStatistic), lower.tail = FALSE)
    }
    colnames(WaldPvalue) <- paste0("WaldPvalue_", modelMatrixNames)
    betaConv <- fit$betaConv
    if (any(!betaConv)) {
        if (!quiet)
            message(paste(sum(!betaConv), "rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest"))
    }
    mleBetas <- if (betaPrior) {
        matrixToList(mleBetaMatrix)
    }
    else {
        NULL
    }
    tDFList <- if (useT)
        list(tDegreesFreedom = df)
    else NULL
    resultsList <- c(matrixToList(betaMatrix), matrixToList(betaSE),
        mleBetas, matrixToList(WaldStatistic), matrixToList(WaldPvalue),
        list(betaConv = betaConv, betaIter = fit$betaIter, deviance = -2 *
            fit$logLike, maxCooks = maxCooks), tDFList)
    WaldResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
    modelMatrixNamesSpaces <- gsub("_", " ", modelMatrixNames)
    lfcType <- if (attr(object, "betaPrior"))
        "MAP"
    else "MLE"
    coefInfo <- paste(paste0("log2 fold change (", lfcType, "):"),
        modelMatrixNamesSpaces)
    seInfo <- paste("standard error:", modelMatrixNamesSpaces)
    mleInfo <- if (betaPrior) {
        gsub("_", " ", colnames(mleBetaMatrix))
    }
    else {
        NULL
    }
    statInfo <- paste("Wald statistic:", modelMatrixNamesSpaces)
    pvalInfo <- paste("Wald test p-value:", modelMatrixNamesSpaces)
    tDFDescription <- if (useT)
        "t degrees of freedom for Wald test"
    else NULL
    mcolsWaldResults <- DataFrame(type = rep("results", ncol(WaldResults)),
        description = c(coefInfo, seInfo, mleInfo, statInfo,
            pvalInfo, "convergence of betas", "iterations for betas",
            "deviance for the fitted model", "maximum Cook's distance for row",
            tDFDescription))
    mcols(WaldResults) <- mcolsWaldResults
    mcols(object) <- cbind(mcols(object), WaldResults)
    return(object)
}
###################################


> DESeq
function (object, test = c("Wald", "LRT"), fitType = c("parametric",
    "local", "mean", "glmGamPoi"), sfType = c("ratio", "poscounts",
    "iterate"), betaPrior, full = design(object), reduced, quiet = FALSE,
    minReplicatesForReplace = 7, modelMatrixType, useT = FALSE,
    minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5, parallel = FALSE,
    BPPARAM = bpparam())
{
    stopifnot(is(object, "DESeqDataSet"))
    test <- match.arg(test, choices = c("Wald", "LRT"))
    fitType <- match.arg(fitType, choices = c("parametric", "local",
        "mean", "glmGamPoi"))
    dispersionEstimator <- if (fitType == "glmGamPoi") {
        "glmGamPoi"
    }
    else {
        "DESeq2"
    }
    if (fitType == "glmGamPoi") {
        minReplicatesForReplace <- Inf
        if (parallel) {
            warning("parallelization of DESeq() is not implemented for fitType='glmGamPoi'")
        }
    }
    sfType <- match.arg(sfType, choices = c("ratio", "poscounts",
        "iterate"))
    stopifnot(is.logical(quiet))
    stopifnot(is.numeric(minReplicatesForReplace))
    stopifnot(is.logical(parallel))
    modelAsFormula <- !is.matrix(full) & is(design(object), "formula")
    if (missing(betaPrior)) {
        betaPrior <- FALSE
    }
    else {
        stopifnot(is.logical(betaPrior))
    }
    object <- sanitizeRowRanges(object)
    if (test == "LRT") {
        if (missing(reduced)) {
            stop("likelihood ratio test requires a 'reduced' design, see ?DESeq")
        }
        if (betaPrior) {
            stop("test='LRT' does not support use of LFC shrinkage, use betaPrior=FALSE")
        }
        if (!missing(modelMatrixType) && modelMatrixType == "expanded") {
            stop("test='LRT' does not support use of expanded model matrix")
        }
        if (is.matrix(full) | is.matrix(reduced)) {
            if (!(is.matrix(full) & is.matrix(reduced))) {
                stop("if one of 'full' and 'reduced' is a matrix, the other must be also a matrix")
            }
        }
        if (modelAsFormula) {
            checkLRT(full, reduced)
        }
        else {
            checkFullRank(full)
            checkFullRank(reduced)
            if (ncol(full) <= ncol(reduced)) {
                stop("the number of columns of 'full' should be more than the number of columns of 'reduced'")
            }
        }
    }
    if (test == "Wald" & !missing(reduced)) {
        stop("'reduced' ignored when test='Wald'")
    }
    if (dispersionEstimator == "glmGamPoi" && test == "Wald") {
        warning("glmGamPoi dispersion estimator should be used in combination with a LRT and not a Wald test.",
            call. = FALSE)
    }
    if (modelAsFormula) {
        designAndArgChecker(object, betaPrior)
        if (design(object) == formula(~1)) {
            warning("the design is ~ 1 (just an intercept). is this intended?")
        }
        if (full != design(object)) {
            stop("'full' specified as formula should equal design(object)")
        }
        modelMatrix <- NULL
    }
    else {
        if (!quiet)
            message("using supplied model matrix")
        if (betaPrior == TRUE) {
            stop("betaPrior=TRUE is not supported for user-provided model matrices")
        }
        checkFullRank(full)
        modelMatrix <- full
    }
    attr(object, "betaPrior") <- betaPrior
    stopifnot(length(parallel) == 1 & is.logical(parallel))
    if (!is.null(sizeFactors(object)) || !is.null(normalizationFactors(object))) {
        if (!quiet) {
            if (!is.null(normalizationFactors(object))) {
                message("using pre-existing normalization factors")
            }
            else {
                message("using pre-existing size factors")
            }
        }
    }
    else {
        if (!quiet)
            message("estimating size factors")
        object <- estimateSizeFactors(object, type = sfType,
            quiet = quiet)
    }
    if (!parallel) {
        if (!quiet)
            message("estimating dispersions")
        object <- estimateDispersions(object, fitType = fitType,
            quiet = quiet, modelMatrix = modelMatrix, minmu = minmu)
        if (!quiet)
            message("fitting model and testing")
        if (test == "Wald") {
            object <- nbinomWaldTest(object, betaPrior = betaPrior,
                quiet = quiet, modelMatrix = modelMatrix, modelMatrixType = modelMatrixType,
                useT = useT, minmu = minmu)
        }
        else if (test == "LRT") {
            object <- nbinomLRT(object, full = full, reduced = reduced,
                quiet = quiet, minmu = minmu, type = dispersionEstimator)
        }
    }
    else if (parallel) {
        if (!missing(modelMatrixType)) {
            if (betaPrior)
                stopifnot(modelMatrixType == "expanded")
        }
        object <- DESeqParallel(object, test = test, fitType = fitType,
            betaPrior = betaPrior, full = full, reduced = reduced,
            quiet = quiet, modelMatrix = modelMatrix, useT = useT,
            minmu = minmu, BPPARAM = BPPARAM)
    }
    sufficientReps <- any(nOrMoreInCell(attr(object, "modelMatrix"),
        minReplicatesForReplace))
    if (sufficientReps) {
        object <- refitWithoutOutliers(object, test = test, betaPrior = betaPrior,
            full = full, reduced = reduced, quiet = quiet, minReplicatesForReplace = minReplicatesForReplace,
            modelMatrix = modelMatrix, modelMatrixType = modelMatrixType)
    }
    metadata(object)[["version"]] <- packageVersion("DESeq2")
    object
