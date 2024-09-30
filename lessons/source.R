> getMethod(estimateDispersions, "DESeqDataSet")
Method Definition:

function (object, ...)
{
    .local <- function (object, fitType = c("parametric", "local",
        "mean", "glmGamPoi"), maxit = 100, useCR = TRUE, weightThreshold = 0.01,
        quiet = FALSE, modelMatrix = NULL, minmu = if (fitType ==
            "glmGamPoi")
            1e-06
        else 0.5)
    {
        if (!.hasSlot(object, "rowRanges"))
            object <- updateObject(object)
        if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
            stop("first call estimateSizeFactors or provide a normalizationFactor matrix before estimateDispersions")
        }
        if (!is.null(sizeFactors(object))) {
            if (!is.numeric(sizeFactors(object))) {
                stop("the sizeFactor column in colData is not numeric.\nthis column could have come in during colData import and should be removed.")
            }
            if (any(is.na(sizeFactors(object)))) {
                stop("the sizeFactor column in colData contains NA.\nthis column could have come in during colData import and should be removed.")
            }
        }
        if (all(MatrixGenerics::rowSums(counts(object) == counts(object)[,
            1]) == ncol(object))) {
            stop("all genes have equal values for all samples. will not be able to perform differential analysis")
        }
        if (!is.null(dispersions(object))) {
            if (!quiet)
                message("found already estimated dispersions, replacing these")
            mcols(object) <- mcols(object)[, !(mcols(mcols(object))$type %in%
                c("intermediate", "results")), drop = FALSE]
        }
        stopifnot(length(maxit) == 1)
        fitType <- match.arg(fitType, choices = c("parametric",
            "local", "mean", "glmGamPoi"))
        dispersionEstimator <- if (fitType == "glmGamPoi") {
            if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
                stop("type='glmGamPoi' requires installing the Bioconductor package 'glmGamPoi'")
            }
            "glmGamPoi"
        }else {
            "DESeq2"
        }
        DESeq2:::checkForExperimentalReplicates(object, modelMatrix)
        if (!quiet)
            message("gene-wise dispersion estimates")
        # genewise dispersion estimates
        object <- DESeq2:::estimateDispersionsGeneEst(object, maxit = maxit,
            useCR = useCR, weightThreshold = weightThreshold,
            quiet = quiet, modelMatrix = modelMatrix, minmu = minmu,
            type = dispersionEstimator)
        if (!quiet)
            message("mean-dispersion relationship")
        # fit mean-dispersion relationship global trend
        object <- DESeq2:::estimateDispersionsFit(object, fitType = fitType,
            quiet = quiet)
        if (!quiet)
            message("final dispersion estimates")
        # Shrinking Gene-Wise Dispersion Estimates (MAP)
        object <- DESeq2:::estimateDispersionsMAP(object, maxit = maxit,
            useCR = useCR, weightThreshold = weightThreshold,
            quiet = quiet, modelMatrix = modelMatrix, type = dispersionEstimator)
        return(object)
    }
    .local(object, ...)
}


> DESeq2:::estimateDispersionsFit
function (object, fitType = c("parametric", "local", "mean",
    "glmGamPoi"), minDisp = 1e-08, quiet = FALSE)
{
    if (is.null(mcols(object)$allZero)) {
        object <- DESeq2:::getBaseMeansAndVariances(object)
    }
    objectNZ <- object[!mcols(object)$allZero, , drop = FALSE]
    useForFit <- mcols(objectNZ)$dispGeneEst > 100 * minDisp
    if (sum(useForFit) == 0) {
        stop("all gene-wise dispersion estimates are within 2 orders of magnitude\n  from the minimum value, and so the standard curve fitting techniques will not work.\n  One can instead use the gene-wise estimates as final estimates:\n  dds <- estimateDispersionsGeneEst(dds)\n  dispersions(dds) <- mcols(dds)$dispGeneEst\n  ...then continue with testing using nbinomWaldTest or nbinomLRT")
    }
    fitType <- match.arg(fitType, choices = c("parametric", "local",
        "mean", "glmGamPoi"))
    stopifnot(length(fitType) == 1)
    stopifnot(length(minDisp) == 1)
    if (fitType == "parametric") {
        trial <- try(dispFunction <- DESeq2:::parametricDispersionFit(mcols(objectNZ)$baseMean[useForFit],
            mcols(objectNZ)$dispGeneEst[useForFit]), silent = TRUE)
        if (inherits(trial, "try-error")) {
            message("-- note: fitType='parametric', but the dispersion trend was not well captured by the\n   function: y = a/x + b, and a local regression fit was automatically substituted.\n   specify fitType='local' or 'mean' to avoid this message next time.")
            fitType <- "local"
        }
    }
    if (fitType == "local") {
        dispFunction <- DESeq2:::localDispersionFit(means = mcols(objectNZ)$baseMean[useForFit],
            disps = mcols(objectNZ)$dispGeneEst[useForFit], minDisp = minDisp)
    }
    if (fitType == "mean") {
        useForMean <- mcols(objectNZ)$dispGeneEst > 10 * minDisp
        meanDisp <- mean(mcols(objectNZ)$dispGeneEst[useForMean],
            na.rm = TRUE, trim = 0.001)
        dispFunction <- function(means) meanDisp
        attr(dispFunction, "mean") <- meanDisp
    }
    if (fitType == "glmGamPoi") {
        if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
            stop("type='glmGamPoi' requires installing the Bioconductor package 'glmGamPoi'")
        }
        base_means <- mcols(objectNZ)$baseMean[useForFit]
        median_fit <- glmGamPoi::loc_median_fit(base_means, mcols(objectNZ)$dispGeneEst[useForFit])
        get_closest_index <- function(x, vec) {
            iv <- findInterval(x, vec)
            dist_left <- x - vec[ifelse(iv == 0, NA, iv)]
            dist_right <- vec[iv + 1] - x
            ifelse(!is.na(dist_left) & (is.na(dist_right) | dist_left <
                dist_right), iv, iv + 1)
        }
        sorted_bm <- sort(base_means)
        ordered_medians <- median_fit[order(base_means)]
        dispFunction <- function(means) {
            indices <- get_closest_index(means, sorted_bm)
            ordered_medians[indices]
        }
    }
    if (!(fitType %in% c("parametric", "local", "mean", "glmGamPoi"))) {
        stop("unknown fitType")
    }
    attr(dispFunction, "fitType") <- fitType
    if (quiet) {
        suppressMessages({
            dispersionFunction(object) <- dispFunction
        })
    }
    else {
        dispersionFunction(object) <- dispFunction
    }
    return(object)
}

DESeq2:::estimateDispersionsGeneEst
function (object, minDisp = 1e-08, kappa_0 = 1, dispTol = 1e-06,
    maxit = 100, useCR = TRUE, weightThreshold = 0.01, quiet = FALSE,
    modelMatrix = NULL, niter = 1, linearMu = NULL, minmu = if (type ==
        "glmGamPoi") 1e-06 else 0.5, alphaInit = NULL, type = c("DESeq2",
        "glmGamPoi"))
{
    type <- match.arg(type, c("DESeq2", "glmGamPoi"))
    if (!is.null(mcols(object)$dispGeneEst)) {
        if (!quiet)
            message("found already estimated gene-wise dispersions, removing these")
        removeCols <- c("dispGeneEst", "dispGeneIter")
        mcols(object) <- mcols(object)[, !names(mcols(object)) %in%
            removeCols, drop = FALSE]
    }
    stopifnot(length(minDisp) == 1)
    stopifnot(length(kappa_0) == 1)
    stopifnot(length(dispTol) == 1)
    stopifnot(length(maxit) == 1)
    if (log(minDisp/10) <= -30) {
        stop("for computational stability, log(minDisp/10) should be above -30")
    }
    object <- DESeq2:::sanitizeRowRanges(object)
    if (is.null(modelMatrix)) {
        modelMatrix <- DESeq2:::getModelMatrix(object)
    }
    DESeq2:::checkFullRank(modelMatrix)
    if (nrow(modelMatrix) == ncol(modelMatrix)) {
        stop("the number of samples and the number of model coefficients are equal,\n  i.e., there are no replicates to estimate the dispersion.\n  use an alternate design formula")
    }
    object <- DESeq2:::getBaseMeansAndVariances(object)
    attr(object, "weightsOK") <- NULL
    wlist <- DESeq2:::getAndCheckWeights(object, modelMatrix, weightThreshold = weightThreshold)
    object <- wlist$object
    weights <- wlist$weights
    weights <- pmax(weights, 1e-06)
    useWeights <- wlist$useWeights
    objectNZ <- object[!mcols(object)$allZero, , drop = FALSE]
    weights <- weights[!mcols(object)$allZero, , drop = FALSE]
    if (is.null(alphaInit)) {
        roughDisp <- DESeq2:::roughDispEstimate(y = counts(objectNZ, normalized = TRUE),
            x = modelMatrix)
        momentsDisp <- DESeq2:::momentsDispEstimate(objectNZ)
        alpha_hat <- pmin(roughDisp, momentsDisp)
    }  else {
        if (length(alphaInit) == 1) {
            alpha_hat <- rep(alphaInit, nrow(objectNZ))
        }else {
            stopifnot(length(alphaInit) == nrow(objectNZ))
            alpha_hat <- alphaInit
        }
    }
    maxDisp <- max(10, ncol(object))
    alpha_hat <- alpha_hat_new <- alpha_init <- pmin(pmax(minDisp,
        alpha_hat), maxDisp)
    stopifnot(length(niter) == 1 & niter > 0)
    if (is.null(linearMu)) {
        modelMatrixGroups <- DESeq2:::modelMatrixGroups(modelMatrix)
        linearMu <- nlevels(modelMatrixGroups) == ncol(modelMatrix)
        if (useWeights) {
            linearMu <- FALSE
        }
    }
    fitidx <- rep(TRUE, nrow(objectNZ))
    mu <- matrix(0, nrow = nrow(objectNZ), ncol = ncol(objectNZ))
    dispIter <- numeric(nrow(objectNZ))
    for (iter in seq_len(niter)) {
        if (!linearMu) {
            fit <- DESeq2:::fitNbinomGLMs(objectNZ[fitidx, , drop = FALSE],
                alpha_hat = alpha_hat[fitidx], modelMatrix = modelMatrix,
                type = type)
            fitMu <- fit$mu
        } else {
            fitMu <- DESeq2:::linearModelMuNormalized(objectNZ[fitidx,
                , drop = FALSE], modelMatrix)
        }
        fitMu[fitMu < minmu] <- minmu
        mu[fitidx, ] <- fitMu
        if (type == "DESeq2") {
            dispRes <- DESeq2:::fitDispWrapper(ySEXP = counts(objectNZ)[fitidx,
                , drop = FALSE], xSEXP = modelMatrix, mu_hatSEXP = fitMu,
                log_alphaSEXP = log(alpha_hat)[fitidx], log_alpha_prior_meanSEXP = log(alpha_hat)[fitidx],
                log_alpha_prior_sigmasqSEXP = 1, min_log_alphaSEXP = log(minDisp/10),
                kappa_0SEXP = kappa_0, tolSEXP = dispTol, maxitSEXP = maxit,
                usePriorSEXP = FALSE, weightsSEXP = weights,
                useWeightsSEXP = useWeights, weightThresholdSEXP = weightThreshold,
                useCRSEXP = useCR)
            dispIter[fitidx] <- dispRes$iter
            alpha_hat_new[fitidx] <- pmin(exp(dispRes$log_alpha),
                maxDisp)
            last_lp <- dispRes$last_lp
            initial_lp <- dispRes$initial_lp
        } else if (type == "glmGamPoi") {
            if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
                stop("type='glmGamPoi' requires installing the Bioconductor package 'glmGamPoi'")
            }
            if (!quiet)
                message("using 'glmGamPoi' as fitType. If used in published research, please cite:\n    Ahlmann-Eltze, C., Huber, W. (2020) glmGamPoi: Fitting Gamma-Poisson\n    Generalized Linear Models on Single Cell Count Data. Bioinformatics.\n    https://doi.org/10.1093/bioinformatics/btaa1009")
            Counts <- counts(objectNZ)[fitidx, , drop = FALSE]
            initial_lp <- vapply(seq_len(nrow(fitMu)), function(idx) {
                sum(dnbinom(x = Counts[idx, ], mu = fitMu[idx,
                  ], size = 1/alpha_hat[fitidx][idx], log = TRUE))
            }, FUN.VALUE = 0)
            dispersion_fits <- glmGamPoi::overdispersion_mle(Counts,
                mean = fitMu, model_matrix = modelMatrix, verbose = !quiet)
            dispIter[fitidx] <- dispersion_fits$iterations
            alpha_hat_new[fitidx] <- pmin(dispersion_fits$estimate,
                maxDisp)
            last_lp <- vapply(seq_len(nrow(fitMu)), function(idx) {
                sum(dnbinom(x = Counts[idx, ], mu = fitMu[idx,
                  ], size = 1/alpha_hat_new[fitidx][idx], log = TRUE))
            }, FUN.VALUE = 0)
        }
        fitidx <- abs(log(alpha_hat_new) - log(alpha_hat)) >
            0.05
        fitidx[is.na(fitidx)] <- FALSE
        alpha_hat <- alpha_hat_new
        if (sum(fitidx) == 0)
            break
    }
    dispGeneEst <- alpha_hat
    if (niter == 1) {
        noIncrease <- last_lp < initial_lp + abs(initial_lp)/1e+06
        dispGeneEst[which(noIncrease)] <- alpha_init[which(noIncrease)]
    }

    # break
    dispGeneEstConv <- dispIter < maxit & !(dispIter == 1)
    refitDisp <- !dispGeneEstConv & dispGeneEst > minDisp * 10
    if (sum(refitDisp) > 0) {
        dispGrid <- DESeq2:::fitDispGridWrapper(y = counts(objectNZ)[refitDisp,
            , drop = FALSE], x = modelMatrix, mu = mu[refitDisp,
            , drop = FALSE], logAlphaPriorMean = rep(0, sum(refitDisp)),
            logAlphaPriorSigmaSq = 1, usePrior = FALSE, weightsSEXP = weights[refitDisp,
                , drop = FALSE], useWeightsSEXP = useWeights,
            weightThresholdSEXP = weightThreshold, useCRSEXP = useCR)
        dispGeneEst[refitDisp] <- dispGrid
    }
    dispGeneEst <- pmin(pmax(dispGeneEst, minDisp), maxDisp)
    dispDataFrame <- DESeq2:::buildDataFrameWithNARows(list(dispGeneEst = dispGeneEst,
        dispGeneIter = dispIter), mcols(object)$allZero)
    mcols(dispDataFrame) <- DataFrame(type = rep("intermediate",
        ncol(dispDataFrame)), description = c("gene-wise estimates of dispersion",
        "number of iterations for gene-wise"))
    mcols(object) <- cbind(mcols(object), dispDataFrame)
    assays(object, withDimnames = FALSE)[["mu"]] <- buildMatrixWithNARows(mu,
        mcols(object)$allZero)
    return(object)
}

DESeq2:::linearModelMuNormalized
function (object, x)
{
    cts <- counts(object)
    norm.cts <- counts(object, normalized = TRUE)
    muhat <- DESeq2:::linearModelMu(norm.cts, x)
    nf <- DESeq2:::getSizeOrNormFactors(object)
    muhat * nf
}

DESeq2:::roughDispEstimate
function (y, x)
{
    mu <- DESeq2:::linearModelMu(y, x)
    mu <- matrix(pmax(1, mu), ncol = ncol(mu))
    m <- nrow(x)
    p <- ncol(x)
    est <- MatrixGenerics::rowSums(((y - mu)^2 - mu)/mu^2)/(m -
        p)
    pmax(est, 0)
}

DESeq2:::linearModelMu
function (y, x)
{
    qrx <- qr(x)
    Q <- qr.Q(qrx)
    Rinv <- solve(qr.R(qrx))
    (y %*% Q) %*% t(x %*% Rinv)
}


# convenience function for building larger matrices
# by filling in NA rows
buildMatrixWithNARows <- function(m, NARows) {
  mFull <- matrix(NA, ncol=ncol(m), nrow=length(NARows))
  mFull[!NARows,] <- m
  mFull
}

# convenience function for building larger matrices
# by filling in 0 rows
buildMatrixWithZeroRows <- function(m, zeroRows) {
  mFull <- matrix(0, ncol=ncol(m), nrow=length(zeroRows))
  mFull[!zeroRows,] <- m
  mFull
}
# convenience function for building results tables
# out of a list and filling in NA rows
buildDataFrameWithNARows <- function(resultsList, NArows) {
  lengths <- sapply(resultsList,length)
  if (!all(lengths == lengths[1])) {
    stop("lengths of vectors in resultsList must be equal")
  }
  if (sum(!NArows) != lengths[1]) {
    stop("number of non-NA rows must be equal to lengths of vectors in resultsList")
  }
  if (sum(NArows) == 0) {
    return(DataFrame(resultsList))
  }
  dfFull <- DataFrame(lapply(resultsList, function(x) vector(mode(x), length(NArows))))
  dfFull[NArows,] <- NA
  dfFull[!NArows,] <- DataFrame(resultsList)
  dfFull
}
