
#' Using the method of Wickramasuriya et al. (2019), this function (based on Hyndman et al.'s hts library) combines the forecasts at all levels of a hierarchical time series and works for degenerate hierarchies.
#'
#' @param fcasts a vector or a matrix (rows = horizon, columns = ts columns) of forecasts
#' @param Smat a structure matrix detailing the hierarchical structure of the hts. Make sure that the order of the rows align with the order of the forecasts.
#' @param residual a matrix of in-sample residuals (columns = ts columns)
#' @param covariance should a shrinkage estimator or the sample estimator be used? alternatively, a custom covariance matrix can be passed (additionally requires the cov.matrix argument)
#' @param nonnegative not implemented yet.
#' @param cov.type specify how the covariance matrix should be computed (default = complete observations). Note that pairwise.complete.obs may not yield a positive definite matrix!
#' @param cov.matrix specify in case a custom covariance matrix should be used
#'
#' @references
#'  [hts: Hierarchical and Grouped Time Series]<https://cran.r-project.org/package=hts>
#'  [Optimal Forecast Reconciliation for Hierarchical and Grouped Time Series Through Trace Minimization]<doi:10.1080/01621459.2018.1448825>
#'  [Degenerate Hierarchical Time Series Reconciliation With The Minimum Trace Algorithm in R]<doi:10.15488/17729>
#'
#'
#' @return reconciled forecasts
#' @export
#' @example inst/examples/MinT_Example.R
MinT <- function (fcasts, Smat, residual, covariance = c("shr", "sam", "custom"),
                  nonnegative = FALSE,#, algorithms = c("chol", "lu", "cg")
                  cov.type = "complete.obs", cov.matrix = NULL)
{
  if(nonnegative) warning("non-negative forecast reconciliation is not yet implemented here.
                          See the hts library for non-degenerate cases instead.")
  if (is.null(Smat)) {
    stop("Please specify the hierarchical or the grouping structure by providing an Smat.", call. = FALSE)
  }

  #alg <- match.arg(algorithms)
  covar <- match.arg(covariance)
  if(covar == "custom" & is.null(cov.matrix)) stop("cov.matrix has to be defined when covar is 'custom'.")
  res <- residual
  #fcasts <- stats::as.ts(fcasts)
  tspx <- stats::tsp(fcasts)
  cnames <- if(is.matrix(fcasts)) colnames(fcasts) else names(fcasts)

  if (missing(residual))
  {
    stop("MinT needs insample residuals.", call. = FALSE)
  }
  if (covar=="sam")
  {
    n <- nrow(res)
    w.1 <- crossprod(res) / n
    if(is.posdef(w.1)==FALSE)
    {
      stop("MinT needs covariance matrix to be positive definite.", call. = FALSE)
    }
  } else if(covar=="shr"){ # shrinkage
    tar <- lowerD(res)
    shrink <- shrink.estim(res, tar, cov.type = cov.type)
    w.1 <- shrink[[1]]
    lambda <- shrink[[2]]
    if (is.posdef(w.1)==FALSE)
    {
      stop("MinT needs covariance matrix to be positive definite.", call. = FALSE)
    }
  } else{ # custom
    w.1 = cov.matrix
  }

  totalts <- nrow(Smat)
  if (!is.matrix(fcasts)) {
    fcasts <- t(fcasts)
  }
  h <- nrow(fcasts)
  if (ncol(fcasts) != totalts) {
    stop("Argument fcasts requires all the forecasts.", call. = FALSE)
  }

  fcasts <- t(fcasts)

  smat <- as.matrix.csr(Smat)
  if (!is.null(w.1)) {
    w.1 <- as.matrix.csr(w.1)
  }
  allf <- CHOL(fcasts = fcasts, S = smat, weights = w.1, allow.changes = FALSE)

  out = t(allf)
  colnames(out) = cnames
  return(out)
}

